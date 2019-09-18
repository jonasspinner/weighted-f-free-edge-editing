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
    // constexpr static bool restrict_solution_by_k = false;
    struct Config {
        bool add_extended_constraints = true;
        bool sparse_constraints = false;

        int all_lazy = 1;
        int num_threads = 4;
        int timelimit = -1;

        int verbosity = 1;
    } config;

    class FSGCallback : public GRBCallback {
    private:
        Graph m_graph;
        std::unique_ptr<FinderI> m_finder;
        const VertexPairMap<GRBVar> &m_vars;
        VertexPairMap<bool> m_edited;

        Config config;

    public:
        FSGCallback(Options::FSG fsg, Graph graph, const VertexPairMap<GRBVar> &vars, Config config_) :
            m_graph(std::move(graph)), m_vars(vars), m_edited(m_graph.size()), config(config_) {
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

                 if (config.sparse_constraints) {
                    VertexPairMap<bool> used_pairs(m_graph.size());
                    m_finder->find([&](const Subgraph& subgraph) {

                        bool not_covered = false;
                        for (VertexPair uv : subgraph.vertexPairs())
                            if (!used_pairs[uv]) {
                                not_covered = true;
                                break;
                            }


                        if (not_covered) {
                            for (VertexPair uv : subgraph.vertexPairs())
                                used_pairs[uv] = true;
                            ++num_added;
                            addSubgraphConstraint(subgraph);
                        }
                        return false;
                    });
                } else {
                    // Find forbidden subgraphs in current solution and add additional constraints.
                    m_finder->find([&](Subgraph &&subgraph) {
                        addSubgraphConstraint(subgraph);
                        ++num_added;
                        return false;
                    });
                }

                if (config.verbosity)
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



public:
    explicit ILPSolver(Options::FSG fsg) : m_fsg(fsg) {}

    Result solve(Instance instance) override {
        auto &[_, graph, costs] = instance;

        try {
            GRBEnv env = GRBEnv();
            GRBModel model = GRBModel(env);

            model.getEnv().set(GRB_IntParam_LogToConsole, (config.verbosity > 0) ? 1 : 0);
            model.set(GRB_IntParam_LazyConstraints, 1);
            model.set(GRB_IntParam_Threads, config.num_threads);
            if (config.timelimit >= 0)
                model.set(GRB_DoubleParam_TimeLimit, config.timelimit);


            // Add model variables: x_uv == 1  <=>  uv is edited
            VertexPairMap<GRBVar> vars(graph.size());
            for (VertexPair uv : graph.vertexPairs())
                vars[uv] = model.addVar(0.0, 1.0, 0, GRB_BINARY);


            // Set objective. When uv is edited the objective is increased by costs[uv]
            GRBLinExpr obj;
            for (VertexPair uv : graph.vertexPairs())
                obj += costs[uv] * vars[uv];

            model.setObjective(obj, GRB_MINIMIZE);

            size_t num_added = addConstraints(model, vars, graph);

            //if (num_added >= 1'000'000)
            //    return Result::Unsolved();

            // Register callback
            FSGCallback cb(m_fsg, graph, vars, config);
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

    static void coverVariables(VertexPairMap<bool> &covered, const Subgraph &subgraph) {
        for (VertexPair uv : subgraph.vertexPairs())
            covered[uv] = true;
    }

    static bool variablesAreCovered(const VertexPairMap<bool> &covered, const Subgraph &subgraph) {
        for (VertexPair uv : subgraph.vertexPairs())
            if (!covered[uv])
                return false;
        return true;
    }

    size_t addConstraints(GRBModel &model, VertexPairMap<GRBVar> &vars, Graph &graph) {
        auto finder = Finder::make(m_fsg, graph);
        size_t num_added = 0;

        if (config.add_extended_constraints) {
            VertexPairMap<bool> covered(graph.size());

            finder->find([&](Subgraph&& subgraph) {
                for (VertexPair uv : subgraph.vertexPairs())
                    covered[uv] = true;

                GRBLinExpr expr;
                for (VertexPair xy : subgraph.vertexPairs())
                    expr += vars[xy];
                auto constr = model.addConstr(expr >= 1);

                constr.set(GRB_IntAttr_Lazy, config.all_lazy);
                ++num_added;
                return false;
            });

            for (VertexPair uv : graph.vertexPairs()) {
                if (covered[uv]) {
                    graph.toggleEdge(uv);

                    finder->find_near(uv, [&](Subgraph&& subgraph) {

                        GRBLinExpr expr;
                        for (VertexPair xy : subgraph.vertexPairs())
                            if (xy == uv) // xy is an edited pair of vertices
                                expr += (1 - vars[xy]);
                            else
                                expr += vars[xy];
                        auto constr = model.addConstr(expr >= 1);

                        constr.set(GRB_IntAttr_Lazy, config.all_lazy);
                        ++num_added;
                        return false;
                    });

                    graph.toggleEdge(uv);
                }
            }
        } else  {
            finder->find([&](const Subgraph& subgraph) {
                GRBLinExpr expr;
                for (VertexPair xy : subgraph.vertexPairs())
                    expr += vars[xy];
                auto constr = model.addConstr(expr >= 1);

                constr.set(GRB_IntAttr_Lazy, config.all_lazy);
                ++num_added;
                return false;
            });
        }

        if (config.verbosity)
            std::cout << "Extended constraints: added " << num_added << " constraints.\n";
        return num_added;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ILPSOLVER_H
