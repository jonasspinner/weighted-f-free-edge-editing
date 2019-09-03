//
// Created by jonas on 03.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LINEARPROGRAMLOWERBOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LINEARPROGRAMLOWERBOUND_H


#include <gurobi_c++.h>


class LinearProgramLowerBound : public LowerBoundI {
    const Graph &m_graph;
    const VertexPairMap<Cost> &m_costs;
    const VertexPairMap<bool> &m_marked;
    std::unique_ptr<GRBEnv> m_env;
    std::unique_ptr<GRBModel> m_model;
    VertexPairMap<GRBVar> m_variables;

    bool verbose = false;

public:
    /**
     * The code is adapted from Michael Hamann.
     */
    LinearProgramLowerBound(const Instance &instance,
                            const VertexPairMap<bool> &forbidden, std::shared_ptr<FinderI> finder_ref) :
            LowerBoundI(std::move(finder_ref)), m_graph(instance.graph),
            m_costs(instance.costs), m_marked(forbidden), m_env(std::make_unique<GRBEnv>()),
            m_variables(m_graph.size()) {}

    /**
     * Initializes the model.
     */
    void initialize() override {
        try {
            m_model = std::make_unique<GRBModel>(*m_env);
            m_model->set(GRB_IntParam_Threads, 1);
            m_model->getEnv().set(GRB_IntParam_OutputFlag, 1);
            m_model->getEnv().set(GRB_IntParam_LogToConsole, verbose ? 1 : 0);
            GRBLinExpr objective = 0;

            for (VertexPair uv : m_graph.vertexPairs()) {
                m_variables[uv] = m_model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
                if (m_graph.has_edge(uv)) {
                    objective += (1 - m_variables[uv]) * m_costs[uv];
                } else {
                    objective += m_variables[uv] * m_costs[uv];
                }
            }

            m_model->setObjective(objective, GRB_MINIMIZE);

            add_forbidden_subgraphs();
        } catch (GRBException &e) {
            std::cout << e.getMessage() << std::endl;
            throw e;
        }

    }

    /**
     * Fixes the vertex pair and adds all nearby forbidden subgraphs as constraints.
     *
     * @param uv
     */
    void after_mark_and_edit(VertexPair uv) override {
        fix_pair(uv, m_graph.has_edge(uv));

        finder->find_near(uv, [&](const Subgraph &subgraph) {
            add_constraint(subgraph);
            return false;
        });
    }

    /**
     * Fixes the vertex pair to the current state.
     *
     * @param uv
     */
    void after_mark(VertexPair uv) override {
        fix_pair(uv, m_graph.has_edge(uv));
    }

    /**
     * Returns a lower bound.
     * The methods builds a new objective function and solves the model.
     *
     * @return
     */
    Cost result(Cost /*k*/) override {
        GRBLinExpr objective = 0;

        for (VertexPair uv : m_graph.vertexPairs()) {
            // TODO: overwrites changes made by after_mark and after_mark_and_edit
            if (m_marked[uv]) {
                fix_pair(uv, m_graph.has_edge(uv));
            } else {
                relax_pair(uv);
            }

            if (m_graph.has_edge(uv)) {
                objective += (1 - m_variables[uv]) * m_costs[uv];
            } else {
                objective += m_variables[uv] * m_costs[uv];
            }
        }

        m_model->setObjective(objective, GRB_MINIMIZE);

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
            if (verbose)
                std::cout << "found_objective: " << found_objective << " rounded result: " << result << std::endl;
            result = std::floor(found_objective);
        }

#ifndef NDEBUG
        std::vector<VertexPair> edits;
        for (VertexPair uv : m_graph.vertexPairs())
            if (m_graph.has_edge(uv) != (m_variables[uv].get(GRB_DoubleAttr_X) >= 0.99))
                edits.push_back(uv);

        Cost sum = 0;
        for (VertexPair uv : edits)
            sum += m_costs[uv];
        if (verbose) {
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
    void add_constraint(const Subgraph &fs) {
        //GRBLinExpr expr = 0;
        //for (VertexPair uv : fs.vertexPairs())
        //    if (m_graph.has_edge(uv))
        //        expr += (1 - variables[uv]);
        //    else
        //        expr += variables[uv];

        GRBLinExpr expr = 3;
        expr -= m_variables[{fs[0], fs[1]}];
        expr -= m_variables[{fs[1], fs[2]}];
        expr -= m_variables[{fs[2], fs[3]}];
        expr += m_variables[{fs[0], fs[2]}];
        expr += m_variables[{fs[1], fs[3]}];

        m_model->addConstr(expr >= 1);
    }

    /**
     * Fixes the variable of uv to be 1.0 if value is true and 0.0 otherwise.
     *
     * @param uv
     * @param exists
     */
    void fix_pair(VertexPair uv, bool exists) {
        double value = exists ? 1.0 : 0.0;
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
     * Adds all found forbidden subgraphs as constraints to the model.
     *
     * @return
     */
    size_t add_forbidden_subgraphs() {
        size_t num_found = 0;

        finder->find([&](const Subgraph &fs) {
            ++num_found;
            add_constraint(fs);
            return false;
        });

        if (verbose)
            std::cout << "added " << num_found << " constraints" << std::endl;
        return num_found;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_LINEARPROGRAMLOWERBOUND_H
