//
// Created by jonas on 19.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ILPSOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ILPSOLVER_H


#include <gurobi_c++.h>

#include <utility>

#include "../finder/finder_utils.h"
#include "Solver.h"
#include "../Instance.h"
#include "../Configuration.h"


class ILPSolver : public Solver {
private:
    // constexpr static bool restrict_solution_by_k = false;

    class FSGCallback : public GRBCallback {
    private:
        Graph m_graph;
        std::unique_ptr<FinderI> m_finder;
        const VertexPairMap<GRBVar> &m_vars;
        VertexPairMap<bool> m_edited;

        const Configuration &m_config;

    public:
        FSGCallback(const Configuration &config, Graph graph, const VertexPairMap<GRBVar> &vars) :
            m_graph(std::move(graph)), m_vars(vars), m_edited(m_graph.size()), m_config(config) {
            m_finder = Finder::make(m_config.forbidden_subgraphs, m_graph);
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

                 if (m_config.sparse_constraints) {
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

                if (m_config.verbosity)
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

    Configuration m_config;
    GRBEnv m_env;

public:
    explicit ILPSolver(Configuration config) : m_config(std::move(config)), m_env() {}

    Result solve(Instance instance) override {
        auto &[graph, costs, _1, _2, _3] = instance;

        try {
            GRBModel model = GRBModel(m_env);

            model.getEnv().set(GRB_IntParam_LogToConsole, (m_config.verbosity > 0) ? 1 : 0);
            model.set(GRB_IntParam_LazyConstraints, 1);
            model.set(GRB_IntParam_Threads, m_config.num_threads);
            if (m_config.timelimit >= 0)
                model.set(GRB_DoubleParam_TimeLimit, m_config.timelimit);


            // Add model variables: x_uv == 1  <=>  uv is edited
            VertexPairMap<GRBVar> vars(graph.size());
            for (VertexPair uv : graph.vertexPairs())
                vars[uv] = model.addVar(0.0, 1.0, 0, GRB_BINARY);


            // Set objective. When uv is edited the objective is increased by costs[uv]
            GRBLinExpr obj;
            for (VertexPair uv : graph.vertexPairs())
                obj += costs[uv] * vars[uv];

            model.setObjective(obj, GRB_MINIMIZE);

            addConstraints(model, vars, graph);

            //if (num_added >= 1'000'000)
            //    return Result::Unsolved();

            // Register callback
            FSGCallback cb(m_config, graph, vars);
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
        auto finder = Finder::make(m_config.forbidden_subgraphs, graph);
        size_t num_added = 0;

        if (m_config.extended_constraints) {
            VertexPairMap<bool> covered(graph.size());

            finder->find([&](Subgraph&& subgraph) {
                for (VertexPair uv : subgraph.vertexPairs())
                    covered[uv] = true;

                GRBLinExpr expr;
                for (VertexPair xy : subgraph.vertexPairs())
                    expr += vars[xy];
                auto constr = model.addConstr(expr >= 1);

                constr.set(GRB_IntAttr_Lazy, m_config.lazy);
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

                        constr.set(GRB_IntAttr_Lazy, m_config.lazy);
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

                constr.set(GRB_IntAttr_Lazy, m_config.lazy);
                ++num_added;
                return false;
            });
        }

        if (m_config.verbosity)
            std::cout << "Extended constraints: added " << num_added << " constraints.\n";

        return num_added;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ILPSOLVER_H
