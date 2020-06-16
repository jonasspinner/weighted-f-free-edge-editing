#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LPRELAXATION_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LPRELAXATION_H


#include <gurobi_c++.h>

#include <utility>

#include "../finder/finder_utils.h"
#include "LowerBoundI.h"
#include "../Instance.h"


namespace lower_bound {

    class LPRelaxation : public LowerBoundI {
        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;
        std::unique_ptr<GRBEnv> m_env;
        std::unique_ptr<GRBModel> m_model;
        VertexPairMap<GRBVar> m_variables;
        Options::FSG m_fsg;
        Cost k_initial;
        bool m_shall_solve;
        std::vector<std::vector<GRBConstr>> m_constraint_stack;

        Configuration m_config;

        constexpr static bool variable_means_edit = false;

        VertexPairMap<bool> m_edited;

        std::shared_ptr<FinderI> finder;
    public:
        /**
         * The code is adapted from Michael Hamann.
         *
         * This lower bound algorithm is based on the following linear program relaxation of the original problem:
         *
         *      \min \sum_{uv} c(uv) \cdot x_{uv}
         *      s.t.   (1 - x_{uv}) + (1 - x_{va}) + (1 - x_{ab}) + x_{ua} + x_{vb} >= 1     \forall (u, v, a, b) \in \mathcal{F}
         *
         *      \forall u, v \in V: x_{uv} = 1 \iff uv \in E'
         */
        LPRelaxation(const Instance &instance, const VertexPairMap<bool> &forbidden,
                     Configuration config, std::shared_ptr<FinderI> finder_ref) :
                m_graph(instance.graph),
                m_costs(instance.costs),
                m_marked(forbidden),
                m_env(std::make_unique<GRBEnv>()),
                m_variables(m_graph.size()),
                m_fsg(finder->forbidden_subgraphs()),
                k_initial(0),
                m_shall_solve(true),
                m_config(std::move(config)),
                m_edited(m_costs.size(), false),
                finder(std::move(finder_ref)) {
            if (finder->forbidden_subgraphs() != Options::FSG::C4P4) {
                throw std::runtime_error("Only C4P4 allowed as forbidden subgraphs.");
            }
        }

        /**
         * Initializes the model.
         */
        void initialize(Cost k) override {
            k_initial = k;
            m_shall_solve = true;

            try {
                m_model = std::make_unique<GRBModel>(*m_env);
                m_model->set(GRB_IntParam_Threads, 1);
                m_model->getEnv().set(GRB_IntParam_OutputFlag, 1);
                m_model->getEnv().set(GRB_IntParam_LogToConsole, (m_config.verbosity > 0) ? 1 : 0);
                if (m_config.timelimit >= 0)
                    m_model->set(GRB_DoubleParam_TimeLimit, m_config.timelimit);

                GRBLinExpr objective = 0;
                for (VertexPair uv : m_graph.vertexPairs()) {
                    m_variables[uv] = m_model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
                    if (variable_means_edit) {
                        objective += m_variables[uv] * m_costs[uv];
                    } else {
                        if (m_graph.hasEdge(uv)) {
                            objective += (1 - m_variables[uv]) * m_costs[uv];
                        } else {
                            objective += m_variables[uv] * m_costs[uv];
                        }
                    }
                }

                m_model->setObjective(objective, GRB_MINIMIZE);

                add_constraints_for_all_forbidden_subgraphs();

                if (variable_means_edit)
                    m_model->addConstr(objective <= k_initial);
            } catch (GRBException &e) {
                std::cerr << e.getMessage() << std::endl;
                std::stringstream ss;
                ss << "GRBException errorCode: " << e.getErrorCode() << ", message: " << e.getMessage();
                throw std::runtime_error(ss.str());
            }

        }

        /**
         * Fixes the vertex pair to the current state.
         *
         * @param uv
         */
        void after_mark(VertexPair uv) override {
            if (!variable_means_edit)
                fix_pair(uv, m_graph.hasEdge(uv));
        }

        void after_edit(VertexPair uv) override {
            if (variable_means_edit) {
                m_edited[uv] = !m_edited[uv];
                fix_pair(uv, m_edited[uv]);
            }
            if (variable_means_edit) {
                assert(m_edited[uv]);
                fix_pair(uv, m_edited[uv]);
            } else {
                fix_pair(uv, m_graph.hasEdge(uv));
            }

            m_constraint_stack.emplace_back();

            if (m_shall_solve) {
                assert(k_initial > 0);
                if (solve() > k_initial) return;
                m_shall_solve = false;
            }

            finder->find_near(uv, m_graph, [&](const Subgraph &subgraph) {
                m_constraint_stack.back().push_back(add_constraint(subgraph));
                if (!m_shall_solve && get_constraint_value(subgraph) < 0.999) {
                    m_shall_solve = true;
                }

                return false;
            });
        }

        void after_unedit(VertexPair uv) override {
            if (variable_means_edit) {
                m_edited[uv] = !m_edited[uv];
                fix_pair(uv, m_edited[uv]);
            }

            fix_pair(uv, m_graph.hasEdge(uv));

            for (auto constr : m_constraint_stack.back()) {
                m_model->remove(constr);
            }

            m_constraint_stack.pop_back();
        }

        void after_unmark(VertexPair uv) override {
            assert(!m_edited[uv]);
            if (variable_means_edit) {
                assert(!m_edited[uv]);
                relax_pair(uv);
            }
            relax_pair(uv);
        }

        /**
         * Returns a lower bound.
         * The methods builds a new objective function and solves the model.
         *
         * @return
         */
        Cost calculate_lower_bound(Cost k) override {
            // /*
            GRBLinExpr objective = 0;

            for (VertexPair uv : m_graph.vertexPairs()) {
                // TODO: overwrites changes made by after_mark and after_mark_and_edit
                if (m_marked[uv]) {
                    if (variable_means_edit)
                        fix_pair(uv, m_edited[uv]);
                    else
                        fix_pair(uv, m_graph.hasEdge(uv));
                } else {
                    relax_pair(uv);
                }

                if (!variable_means_edit) {
                    if (m_graph.hasEdge(uv)) {
                        objective += (1 - m_variables[uv]) * m_costs[uv];
                    } else {
                        objective += m_variables[uv] * m_costs[uv];
                    }
                }
            }

            if (!variable_means_edit)
                m_model->setObjective(objective, GRB_MINIMIZE);

            if (variable_means_edit)
                return solve() - k;
            else
                return solve();
            //*/

            if (m_shall_solve) {
                auto result = solve();
                if (result > k_initial) {
                    return result;
                }
                m_shall_solve = false;
            }
            return 0;
        }

    public:
        [[nodiscard]] auto variable_edge_value(VertexPair uv) const {
            return m_variables[uv].get(GRB_DoubleAttr_X);
        }

        [[nodiscard]] auto variable_edited_value(VertexPair uv) const {
            if (m_graph.hasEdge(uv)) {
                return 1 - m_variables[uv].get(GRB_DoubleAttr_X);
            } else {
                return m_variables[uv].get(GRB_DoubleAttr_X);
            }
        }

    private:
        /**
         * Solves the current model and returns the objective value.
         *
         * @return
         */
        Cost solve() {
            m_model->optimize();
            if (m_model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
                return std::numeric_limits<Cost>::max();
            }
            assert(m_model->get(GRB_IntAttr_Status) == GRB_OPTIMAL);

            double found_objective = m_model->get(GRB_DoubleAttr_ObjVal);
            Cost result = std::ceil(found_objective);

            if (result - found_objective > 0.99) {
                if (m_config.verbosity)
                    std::cout << "found_objective: " << found_objective << " rounded result: " << result << std::endl;
                result = std::floor(found_objective);
            }

#ifndef NDEBUG
            std::vector<VertexPair> edits;
            for (VertexPair uv : m_graph.vertexPairs())
                if (variable_means_edit) {
                    if (m_variables[uv].get(GRB_DoubleAttr_X) >= 0.99)
                        edits.push_back(uv);
                } else {
                    if (m_graph.hasEdge(uv) != (m_variables[uv].get(GRB_DoubleAttr_X) >= 0.99))
                        edits.push_back(uv);
                }

            Cost sum = 0;
            for (VertexPair uv : edits)
                sum += m_costs[uv];
            if (m_config.verbosity) {
                std::cout << "lower bound ";
                for (VertexPair uv : edits)
                    std::cout << uv << " " << m_costs[uv] << "  ";
                std::cout << sum << " " << result << "\n";
            }
            // assert(sum == result);
#endif
            return result;
        }

        /**
         * Adds a subgraph as a contraint. Assumes that the subgraph is either a P_4 or a C_4.
         * @param fs
         */
        GRBConstr add_constraint(const Subgraph &subgraph) {
            /*
            Vertex u = subgraph[0], v = subgraph[1], a = subgraph[2], b = subgraph[3];
            auto x = [&](VertexPair e) { return m_variables[e]; };

            GRBLinExpr expr = 3.0 - x({u, v}) - x({v, a}) - x({a, b}) + x({u, a}) + x({v, b});

            return m_model->addConstr(expr >= 1);
             */
            return m_model->addConstr(alpha(subgraph) >= 1);
        }

        /**
         * Fixes the variable of uv to be 1.0 if value is true and 0.0 otherwise.
         *
         * @param uv
         * @param exists
         */
        void fix_pair(VertexPair uv, bool exists) {
            const double epsilon = 0.001;
            double value = exists ? 1.0 : 0.0;

            if (!m_shall_solve) {
                if (exists) {
                    m_shall_solve = m_variables[uv].get(GRB_DoubleAttr_X) < 1 - epsilon;
                } else {
                    m_shall_solve = m_variables[uv].get(GRB_DoubleAttr_X) > epsilon;
                }
            }

            m_variables[uv].set(GRB_DoubleAttr_UB, value);
            m_variables[uv].set(GRB_DoubleAttr_LB, value);
        }

        /**
         * Relaxes the restrictions on the variable of uv. The bounds are set to [0.0, 1.0].
         *
         * @param uv
         */
        void relax_pair(VertexPair uv) {
            m_variables[uv].set(GRB_DoubleAttr_LB, 0.0);
            m_variables[uv].set(GRB_DoubleAttr_UB, 1.0);
        }

        /**
         * Adds all found forbidden subgraphs currently in the graph as constraints to the model.
         *
         * @return
         */
        size_t add_constraints_for_all_forbidden_subgraphs() {
            size_t num_found = 0;

            finder->find_with_duplicates(m_graph, [&](const Subgraph &subgraph) {
                ++num_found;
                add_constraint(subgraph);
                return false;
            });

            if (m_config.verbosity)
                std::cout << "added " << num_found << " constraints" << std::endl;
            return num_found;
        }

        /**
         * The linear program has constraints of the form
         *
         *      $\alpha \cdot x \geq 1$.
         *
         *  This function calculates $\alpha \cdot x$ for a given $x$.
         */
        double get_constraint_value(const Subgraph &subgraph) {
            assertC4orP4(m_graph, subgraph);
            /*
            auto x = [&](VertexPair e) {
                if (m_edited[e]) {
                    return m_graph.hasEdge(e) ? 1.0 : 0.0;
                }
                return m_variables[e].get(GRB_DoubleAttr_X);
            };

            // v---a
            // |   |
            // u-?-b
            Vertex u = subgraph[0], v = subgraph[1], a = subgraph[2], b = subgraph[3];
            return 3.0 - x({u, v}) - x({v, a}) - x({a, b}) + x({u, a}) + x({v, b});
             */
            return alpha(subgraph, true).getValue();
        }


        GRBLinExpr alpha(const Subgraph &subgraph, bool fixed_var_is_constant = false) {
            assertC4orP4(m_graph, subgraph);
            Vertex u = subgraph[0], v = subgraph[1], a = subgraph[2], b = subgraph[3];

            auto x = [&](VertexPair e) -> GRBLinExpr {
                if (fixed_var_is_constant && m_edited[e]) {
                    return m_graph.hasEdge(e) ? 1.0 : 0.0;
                } else {
                    return m_variables[e];
                }
            };

            return 3.0 - x({u, v}) - x({v, a}) - x({a, b}) + x({u, a}) + x({v, b});
        }

        static void assertC4orP4(const Graph &graph, const Subgraph &subgraph) {
#ifndef NDEBUG
            assert(subgraph.size() == 4);
            Vertex u = subgraph[0], v = subgraph[1], a = subgraph[2], b = subgraph[3];
            assert(graph.hasEdge({u, v}));
            assert(graph.hasEdge({v, a}));
            assert(graph.hasEdge({a, b}));
            assert(!graph.hasEdge({u, a}));
            assert(!graph.hasEdge({v, b}));
#endif
        }
    };

}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_LPRELAXATION_H
