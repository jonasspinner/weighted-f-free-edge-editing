#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ILPSOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ILPSOLVER_H


#include <gurobi_c++.h>

#include <utility>

#include "Solver.h"
#include "../Instance.h"
#include "../Configuration.h"


class ILPSolver : public Solver {
private:
    // constexpr static bool restrict_solution_by_k = false;

    template<Options::FSG FSG>
    class FSGCallback : public GRBCallback {
    private:
        using Subgraph = SubgraphT<Options::FSG::C4P4>;

        Graph m_graph;
        Subgraph::Finder m_finder;
        const VertexPairMap<GRBVar> *m_vars;
        VertexPairMap<bool> m_edited;

        int m_verbosity;
        Options::ILPConstraintGeneration m_constraints;

    public:
        FSGCallback(const Configuration &config, Graph graph, const VertexPairMap<GRBVar> &vars) :
                m_graph(std::move(graph)), m_vars(std::addressof(vars)), m_edited(m_graph.size()),
                m_verbosity(config.verbosity), m_constraints(config.constraints) {
            static_assert(FSG == Options::FSG::C4P4, "Only C4P4 currently supported.");
            if (config.forbidden_subgraphs != Options::FSG::C4P4) {
                throw std::runtime_error("Only C4P4 currently supported.");
            }
        };

    protected:
        void callback() override {
            if (where == GRB_CB_MIPSOL) {

                std::vector<VertexPair> edits;

                // List edits of current solution. Apply them to graph.
                for (VertexPair uv : m_graph.vertex_pairs())
                    if (getSolution((*m_vars)[uv]) > 0.5) {
                        edits.push_back(uv);
                        m_graph.toggle_edge(uv);
                        m_edited[uv] = true;
                    } else {
                        m_edited[uv] = false;
                    }


                size_t num_added = 0;

                switch (m_constraints) {
                    case Options::ILPConstraintGeneration::AtMostOnePerVertexPair: {
                        VertexPairMap<bool> used_pairs(m_graph.size());
                        m_finder.find(m_graph, [&](const Subgraph &subgraph) {

                            bool not_covered = false;
                            for (VertexPair uv : subgraph.vertex_pairs())
                                if (!used_pairs[uv]) {
                                    not_covered = true;
                                    break;
                                }


                            if (not_covered) {
                                for (VertexPair uv : subgraph.vertex_pairs())
                                    used_pairs[uv] = true;
                                ++num_added;
                                addSubgraphConstraint(subgraph);
                            }
                            return subgraph_iterators::IterationControl::Continue;
                        });
                    }
                        break;
                    case Options::ILPConstraintGeneration::One: {
                        m_finder.find(m_graph, [&](const Subgraph &subgraph) {
                            addSubgraphConstraint(subgraph);
                            ++num_added;
                            return subgraph_iterators::IterationControl::Break;
                        });
                    }
                        break;
                    case Options::ILPConstraintGeneration::All: {
                        // Find forbidden subgraphs in current solution and add additional constraints.
                        m_finder.find(m_graph, [&](const Subgraph &subgraph) {
                            addSubgraphConstraint(subgraph);
                            ++num_added;
                            return subgraph_iterators::IterationControl::Continue;
                        });
                    }
                        break;
                    default:
                        throw std::runtime_error("invalid constraint generation.");
                }

                if (m_verbosity)
                    std::cout << "Callback constraints: added " << num_added << " constraints.\n";

                // Undo edits.
                for (VertexPair uv : edits)
                    m_graph.toggle_edge(uv);

            }
        }

    private:
        void addSubgraphConstraint(const Subgraph &subgraph) {
            GRBLinExpr pairs;
            for (VertexPair uv : subgraph.vertex_pairs()) {
                if (m_edited[uv]) {
                    pairs += (1 - (*m_vars)[uv]);
                } else {
                    pairs += (*m_vars)[uv];
                }
            }

            addLazy(pairs >= 1);
        }
    };

    Configuration m_config;
    int m_verbosity;
    int m_timelimit;
    int m_num_threads;
    GRBEnv m_env;

public:
    explicit ILPSolver(Configuration config) : m_config(std::move(config)), m_verbosity(m_config.verbosity),
                                               m_timelimit(m_config.timelimit), m_num_threads(m_config.num_threads),
                                               m_env() {
        if (m_config.forbidden_subgraphs != Options::FSG::C4P4) {
            throw std::runtime_error("Only C4P4 currently supported.");
        }
    }

    Result solve(const Instance &instance) override {
        auto &[graph, costs, _1, _2, _3] = instance;

        try {
            GRBModel model = GRBModel(m_env);

            model.getEnv().set(GRB_IntParam_LogToConsole, (m_verbosity > 0) ? 1 : 0);
            model.set(GRB_IntParam_Threads, m_num_threads);
            model.set(GRB_IntParam_LazyConstraints, 1);
            if (m_timelimit >= 0)
                model.set(GRB_DoubleParam_TimeLimit, m_timelimit);


            // Add model variables: x_uv == 1  <=>  uv is edited
            VertexPairMap<GRBVar> vars(graph.size());
            for (VertexPair uv : graph.vertex_pairs())
                vars[uv] = model.addVar(0.0, 1.0, 0, GRB_BINARY);


            // Set objective. When uv is edited the objective is increased by costs[uv]
            GRBLinExpr obj;
            for (VertexPair uv : graph.vertex_pairs())
                obj += costs[uv] * vars[uv];

            model.setObjective(obj, GRB_MINIMIZE);

            addConstraints<Options::FSG::C4P4>(model, vars, graph);

            // Register callback
            FSGCallback<Options::FSG::C4P4> cb(m_config, graph.copy(), vars);
            model.setCallback(&cb);


            // Optimize model
            model.optimize();

            if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
                return Result::Timeout();


            // Read solution
            std::vector<VertexPair> edits;
            for (VertexPair uv : graph.vertex_pairs())
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

    template<Options::FSG FSG>
    size_t addConstraints(GRBModel &model, const VertexPairMap<GRBVar> &vars, const Graph &graph) {
        using Subgraph = SubgraphT<FSG>;
        typename Subgraph::Finder finder;
        size_t num_added = 0;

        finder.find(graph, [&](const Subgraph &subgraph) {
            GRBLinExpr expr;
            for (VertexPair xy : subgraph.vertex_pairs())
                expr += vars[xy];
            model.addConstr(expr >= 1);

            ++num_added;
            return subgraph_iterators::IterationControl::Continue;
        });

        return num_added;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ILPSOLVER_H
