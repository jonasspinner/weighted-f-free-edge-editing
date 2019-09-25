//
// Created by jonas on 03.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LPRELAXATIONLOWERBOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LPRELAXATIONLOWERBOUND_H


#include <gurobi_c++.h>

#include "../finder/finder_utils.h"


class LPRelaxationLowerBound : public LowerBoundI {
    const Graph &m_graph;
    const VertexPairMap<Cost> &m_costs;
    const VertexPairMap<bool> &m_marked;
    std::unique_ptr<GRBEnv> m_env;
    std::unique_ptr<GRBModel> m_model;
    VertexPairMap<GRBVar> m_variables;
    Options::FSG m_fsg;
    Cost k_initial;

    int verbosity = 0;

    constexpr static bool variable_means_edit = false;

    VertexPairMap<bool> m_edited;

public:
    /**
     * The code is adapted from Michael Hamann.
     */
    LPRelaxationLowerBound(const Instance &instance,
                           const VertexPairMap<bool> &forbidden, std::shared_ptr<FinderI> finder_ref) :
            LowerBoundI(std::move(finder_ref)), m_graph(instance.graph),
            m_costs(instance.costs), m_marked(forbidden), m_env(std::make_unique<GRBEnv>()),
            m_variables(m_graph.size()), m_fsg(finder->forbidden_subgraphs()), k_initial(0), m_edited(m_costs.size(), false) {}

    /**
     * Initializes the model.
     */
    void initialize(Cost k) override {
        k_initial = k;
        try {
            m_model = std::make_unique<GRBModel>(*m_env);
            m_model->set(GRB_IntParam_Threads, 1);
            m_model->getEnv().set(GRB_IntParam_OutputFlag, 1);
            m_model->getEnv().set(GRB_IntParam_LogToConsole, (verbosity > 0) ? 1 : 0);

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

            addForbiddenSubgraphs();

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
     * Fixes the vertex pair and adds all nearby forbidden subgraphs as constraints.
     *
     * @param uv
     */
    void after_mark_and_edit(VertexPair uv) override {
        if (!re_structure) {
            if (variable_means_edit) {
                assert(m_edited[uv]);
                fixPair(uv, m_edited[uv]);
            }
            else
                fixPair(uv, m_graph.hasEdge(uv));

            finder->find_near(uv, [&](const Subgraph &subgraph) {
                addConstraint(subgraph);
                return false;
            });
        }
    }

    /**
     * Fixes the vertex pair to the current state.
     *
     * @param uv
     */
    void after_mark(VertexPair uv) override {
        if (!variable_means_edit)
            fixPair(uv, m_graph.hasEdge(uv));
    }

    void after_edit(VertexPair uv) override {
        if (re_structure) {
            if (variable_means_edit) {
                m_edited[uv] = !m_edited[uv];
                fixPair(uv, m_edited[uv]);
            }

            if (variable_means_edit) {
                assert(m_edited[uv]);
                fixPair(uv, m_edited[uv]);
            }
            else
                fixPair(uv, m_graph.hasEdge(uv));

            finder->find_near(uv, [&](const Subgraph &subgraph) {
                addConstraint(subgraph);
                return false;
            });
        } else {
            // toggle edited status
            if (variable_means_edit) {
                m_edited[uv] = !m_edited[uv];
                fixPair(uv, m_edited[uv]);
            }
        }
    }

    void after_unedit(VertexPair uv) override {
        if constexpr (re_structure) {
            if (variable_means_edit) {
                m_edited[uv] = !m_edited[uv];
                fixPair(uv, m_edited[uv]);
            }
        } else {
            after_edit(uv);
        }
    }

    void after_unmark(VertexPair uv) override {
        if (variable_means_edit) {
            assert(!m_edited[uv]);
            relaxPair(uv);
        }
    }

    /**
     * Returns a lower bound.
     * The methods builds a new objective function and solves the model.
     *
     * @return
     */
    Cost calculate_lower_bound(Cost k) override {
        GRBLinExpr objective = 0;

        for (VertexPair uv : m_graph.vertexPairs()) {
            // TODO: overwrites changes made by after_mark and after_mark_and_edit
            if (m_marked[uv]) {
                if (variable_means_edit)
                    fixPair(uv, m_edited[uv]);
                else
                    fixPair(uv, m_graph.hasEdge(uv));
            } else {
                relaxPair(uv);
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
            if (verbosity)
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
        if (verbosity) {
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
    void addConstraint(const Subgraph &fs) {
        GRBLinExpr expr;
        if (variable_means_edit) {
            for (VertexPair uv : fs.vertexPairs())
                if (m_marked[uv])
                    expr += (1 - m_variables[uv]);
                else
                    expr += m_variables[uv];
        } else {
            expr += 3;
            expr -= m_variables[{fs[0], fs[1]}];
            expr -= m_variables[{fs[1], fs[2]}];
            expr -= m_variables[{fs[2], fs[3]}];
            expr += m_variables[{fs[0], fs[2]}];
            expr += m_variables[{fs[1], fs[3]}];
        }

        m_model->addConstr(expr >= 1);
    }

    /**
     * Fixes the variable of uv to be 1.0 if value is true and 0.0 otherwise.
     *
     * @param uv
     * @param exists
     */
    void fixPair(VertexPair uv, bool exists) {
        double value = exists ? 1.0 : 0.0;
        m_variables[uv].set(GRB_DoubleAttr_UB, value);
        m_variables[uv].set(GRB_DoubleAttr_LB, value);
    }

    /**
     * Relaxes the restrictions on the variable of uv. The bounds are set to [0.0, 1.0].
     *
     * @param uv
     */
    void relaxPair(VertexPair uv) {
        m_variables[uv].set(GRB_DoubleAttr_LB, 0.0);
        m_variables[uv].set(GRB_DoubleAttr_UB, 1.0);
    }

    /**
     * Adds all found forbidden subgraphs currently in the graph as constraints to the model.
     *
     * @return
     */
    size_t addForbiddenSubgraphs() {
        size_t num_found = 0;

        finder->find([&](const Subgraph &fs) {
            ++num_found;
            addConstraint(fs);
            return false;
        });

        if (verbosity)
            std::cout << "added " << num_found << " constraints" << std::endl;
        return num_found;
    }

    void addDistanceOneGraphForbiddenSubgraphs() {
        Graph G(m_graph);
        auto G_finder = Finder::make(m_fsg, G);
        for (VertexPair uv : G.vertexPairs()) {
            G.toggleEdge(uv);

            G_finder->find_near(uv, [&](Subgraph &&subgraph) {
                addConstraint(subgraph);
                return false;
            });

            G.toggleEdge(uv);
        }
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_LPRELAXATIONLOWERBOUND_H
