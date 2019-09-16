//
// Created by jonas on 19.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ILPSOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ILPSOLVER_H


#include <gurobi_c++.h>

#include "../finder/finder_utils.h"
#include "Solver.h"
#include "../Instance.h"
#include "../Configuration.h"


class ILPSolver : public Solver {
private:
    class FSGCallback : public GRBCallback {
    private:
        Graph m_graph;
        std::unique_ptr<FinderI> m_finder;
        const VertexPairMap<GRBVar> &m_vars;
        VertexPairMap<bool> m_edited;

        int verbosity;

    public:
        FSGCallback(Options::FSG fsg, Graph graph, const VertexPairMap<GRBVar> &vars, int verbosity_) :
            m_graph(std::move(graph)), m_vars(vars), m_edited(m_graph.size()), verbosity(verbosity_) {
            m_finder = Finder::make(fsg, m_graph);
        };

    protected:
        void callback() override {
            if (where == GRB_CB_MIPSOL) {

                std::vector<VertexPair> edits;

                // List edits of current solution. Apply them to graph.
                for (VertexPair uv : m_graph.vertexPairs())
                    if (getSolution(m_vars[uv]) > 0.5) {
                        edits.push_back(uv);
                        m_graph.toggleEdge(uv);
                        m_edited[uv] = true;
                    } else {
                        m_edited[uv] = false;
                    }


                size_t num_added = 0;
                // Find forbidden subgraphs in current solution and add additional constraints.
                m_finder->find([&](Subgraph &&subgraph) {
                    addSubgraphConstraint(subgraph);
                    ++num_added;
                    return false;
                });

                if (verbosity)
                    std::cout << "Callback constraints: added " << num_added << " constraints.\n";

                // Undo edits.
                for (VertexPair uv : edits)
                    m_graph.toggleEdge(uv);

            }
        }

    private:
        void addSubgraphConstraint(const Subgraph &subgraph) {
            GRBLinExpr pairs;
            for (VertexPair uv : subgraph.vertexPairs()) {
                if (m_edited[uv]) {
                    pairs += (1 - m_vars[uv]);
                } else {
                    pairs += m_vars[uv];
                }
            }

            addLazy(pairs >= 1);
        }
    };

    Options::FSG m_fsg;


    constexpr static bool restrict_solution_by_k = false;
    constexpr static bool add_initial_constraints = true;
    constexpr static bool add_extended_constraints = true;

    int all_lazy = 1;
    int num_threads = 4;

    int verbosity = 1;

public:
    explicit ILPSolver(Options::FSG fsg) : m_fsg(fsg) {}

    Result solve(Instance instance) override {
        auto &[_, graph, costs] = instance;

        try {
            GRBEnv env = GRBEnv();
            GRBModel model = GRBModel(env);

            model.getEnv().set(GRB_IntParam_LogToConsole, (verbosity > 0) ? 1 : 0);
            model.set(GRB_IntParam_LazyConstraints, 1);
            model.set(GRB_IntParam_Threads, num_threads);
            model.set(GRB_DoubleParam_TimeLimit, 100);


            // Add model variables: x_uv == 1  <=>  uv is edited
            VertexPairMap<GRBVar> vars(graph.size());
            for (VertexPair uv : graph.vertexPairs())
                vars[uv] = model.addVar(0.0, 1.0, 0, GRB_BINARY);


            // Set objective. When uv is edited the objective is increased by costs[uv]
            GRBLinExpr obj;
            for (VertexPair uv : graph.vertexPairs())
                obj += costs[uv] * vars[uv];

            model.setObjective(obj, GRB_MINIMIZE);


            size_t num_added = 0;
            if (add_initial_constraints)
                num_added += addInitialConstraints(model, vars, graph);
            if (add_extended_constraints)
                num_added += addExtendedConstraints(model, vars, graph);

            if (num_added >= 1'000'000)
                return {};

            // Register callback
            FSGCallback cb(m_fsg, graph, vars, verbosity);
            model.setCallback(&cb);


            // Optimize model
            model.optimize();

            if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
                return Result::Timeout();


            // Read solution
            std::vector<VertexPair> edits;
            for (VertexPair uv : graph.vertexPairs())
                if (is_set(vars, uv))
                    edits.push_back(uv);

            return Result::Solved({Solution(instance, edits)});

        } catch (GRBException &e) {
            std::cerr << "Error " << e.getErrorCode() << " " << e.getMessage();
            abort();
        }
    }

private:
    static bool is_set(const VertexPairMap<GRBVar> &vars, VertexPair uv) {
        double value = vars[uv].get(GRB_DoubleAttr_X);
        assert(0.01 < value || value < 0.99);
        return value > 0.5;
    }

    size_t addInitialConstraints(GRBModel &model, VertexPairMap<GRBVar> &vars, Graph &graph) {
        auto finder = Finder::make(m_fsg, graph);
        size_t num_added = 0;

        finder->find([&](Subgraph&& subgraph) {
            GRBLinExpr expr;
            for (VertexPair xy : subgraph.vertexPairs())
                expr += vars[xy];
            auto constr = model.addConstr(expr >= 1);
            constr.set(GRB_IntAttr_Lazy, all_lazy);
            ++num_added;
            return false;
        });
        if (verbosity)
            std::cout << "Initial constraints: added " << num_added << " constraints.\n";
        return num_added;
    }

    size_t addExtendedConstraints(GRBModel &model, VertexPairMap<GRBVar> &vars, Graph &graph) {
        auto finder = Finder::make(m_fsg, graph);
        size_t num_added = 0;

        for (VertexPair uv : graph.vertexPairs()) {
            graph.toggleEdge(uv);

            finder->find_near(uv, [&](Subgraph&& subgraph) {
                GRBLinExpr expr;
                for (VertexPair xy : subgraph.vertexPairs())
                    if (xy == uv)
                        expr += (1 - vars[xy]);
                    else
                        expr += vars[xy];
                auto constr = model.addConstr(expr >= 1);
                constr.set(GRB_IntAttr_Lazy, all_lazy);
                ++num_added;
                return false;
            });

            graph.toggleEdge(uv);
        }
        if (verbosity)
            std::cout << "Extended constraints: added " << num_added << " constraints.\n";
        return num_added;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ILPSOLVER_H
