//
// Created by jonas on 14.01.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GUROBILOWERBOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GUROBILOWERBOUND_H


#include <gurobi_c++.h>

#include "../LowerBoundI.h"
#include "../../Instance.h"


namespace LowerBound {
    class GurobiLowerBound : public LowerBoundI {
    private:
        std::unique_ptr<GRBEnv> env;
        std::unique_ptr<GRBModel> model;
        VertexPairMap<GRBVar> variables;
        Cost initial_k;
        Cost objective_offset;
        bool shall_solve;
        std::vector<std::vector<GRBConstr>> constraint_stack;


        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;

        size_t solve() {
            model->optimize();
            if (model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
                return std::numeric_limits<size_t>::max();
            }

            if (model->get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
                std::cerr << model->get(GRB_IntAttr_Status) << std::endl;
                throw std::runtime_error("model not optimal after optimization");
            }

            double found_objective = model->get(GRB_DoubleAttr_ObjVal);
            size_t result = std::ceil(found_objective);
            if (result - found_objective > 0.99) {
                //std::cerr << "found_objective: " << found_objective << " rounded result: " << result << std::endl;
                result = std::floor(found_objective);
            }

            return result;
        }

        static double get_constraint_value(const Subgraph &subgraph, const Graph &graph, const VertexPairMap<bool> &marked, const VertexPairMap<GRBVar> &variables) {
            auto get = [&](Vertex i, Vertex j) {
                if (marked[{i, j}]) {
                    return graph.hasEdge({i, j}) ? 1.0 : 0.0;
                }
                return variables[{subgraph[i], subgraph[j]}].get(GRB_DoubleAttr_X);
            };

            return 3.0 - get(0, 1) - get(1, 2) - get(2, 3) + get(0, 2) + get(1, 3);
        }

        static GRBConstr add_constraint(GRBModel &model, const VertexPairMap<GRBVar> &variables, const Subgraph &subgraph) {
            Vertex u = subgraph[0], v = subgraph[1], w = subgraph[2], x = subgraph[3];
            GRBLinExpr expr = 3;
            expr -= variables[{u, v}];
            expr -= variables[{v, w}];
            expr -= variables[{w, x}];
            expr += variables[{u, w}];
            expr += variables[{v, x}];

            return model.addConstr(expr, GRB_GREATER_EQUAL, 1);
        }

        void fix_pair(VertexPair uv, bool exists) {
            GRBVar var(variables[uv]);
            if (exists) {
                if (!shall_solve) {
                    // round to prevent too frequent solving, not solving is not a problem.
                    shall_solve = var.get(GRB_DoubleAttr_X) < 0.999;
                }
                var.set(GRB_DoubleAttr_UB, 1.0);
                var.set(GRB_DoubleAttr_LB, 1.0);
            } else {
                if (!shall_solve) {
                    // round to prevent too frequent solving, not solving is not a problem.
                    shall_solve = var.get(GRB_DoubleAttr_X) > 0.001;
                }
                var.set(GRB_DoubleAttr_LB, 0.0);
                var.set(GRB_DoubleAttr_UB, 0.0);
            }
        }

        void relax_pair(VertexPair uv) {
            // No need to set shall_solve here as either we pruned, then shall_solve is set anyway, or we cannot prune
            GRBVar var(variables[uv]);
            var.set(GRB_DoubleAttr_LB, 0.0);
            var.set(GRB_DoubleAttr_UB, 1.0);
        }

        static size_t add_forbidden_subgraphs(GRBModel &model, const VertexPairMap<GRBVar> &variables, FinderI &finder, const Graph &graph) {
            size_t num_found = 0;

            finder.find(graph, [&](const Subgraph &fs) {
                ++num_found;
                add_constraint(model, variables, fs);
                return false;
            });

            std::cerr << "added " << num_found << " constraints" << std::endl;
            return num_found;
        }

    public:
        GurobiLowerBound(const Instance &instance, const VertexPairMap<bool> &marked,
                         std::shared_ptr<FinderI> finder_ref)
                : LowerBoundI(std::move(finder_ref)),
                  env(std::make_unique<GRBEnv>()),
                  variables(instance.graph.size()),
                  initial_k(0), shall_solve(true),
                  m_graph(instance.graph), m_costs(instance.costs),
                  m_marked(marked) {}

        void initialize(Cost k) override {
            set_initial_k(k);

            try {
                model = std::make_unique<GRBModel>(*env);
                model->set(GRB_IntParam_Threads, 1);
                model->getEnv().set(GRB_IntParam_OutputFlag, 1);
                model->getEnv().set(GRB_IntParam_LogToConsole, 0);
                GRBLinExpr objective = 0;
                objective_offset = 0;

                for (VertexPair uv : m_graph.vertexPairs()) {
                    variables[uv] = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
                    if (m_graph.hasEdge(uv)) {
                        objective -= variables[uv] * m_costs[uv];
                        ++objective_offset;
                    } else {
                        objective += variables[uv] * m_costs[uv];
                    }
                }

                objective += objective_offset;
                model->setObjective(objective, GRB_MINIMIZE);

                add_forbidden_subgraphs(*model, variables, *finder, m_graph);

                shall_solve = true;
            } catch (GRBException &e) {
                std::cerr << e.getMessage() << std::endl;
                throw e;
            }
        }

        void set_initial_k(Cost k) {
            initial_k = k;
            shall_solve = true;
        }

        void after_mark_and_edit(VertexPair uv) override {
            try {
                fix_pair(uv, m_graph.hasEdge(uv));

                constraint_stack.emplace_back();

                if (shall_solve) {
                    assert(initial_k > 0);
                    if (solve() > initial_k) return;
                    shall_solve = false;
                }

                finder->find_near(uv, [&](const Subgraph &subgraph) {
                    constraint_stack.back().push_back(add_constraint(*model, variables, subgraph));
                    if (!shall_solve && get_constraint_value(subgraph, m_graph, m_marked, variables) < 0.999) {
                        shall_solve = true;
                    }

                    return false;
                });
            } catch (GRBException &e) {
                std::cerr << e.getMessage() << std::endl;
                throw e;
            }
        }

        void after_undo_edit(VertexPair uv) override {
            try {
                fix_pair(uv, m_graph.hasEdge(uv));

                for (GRBConstr cstr : constraint_stack.back()) {
                    model->remove(cstr);
                }

                constraint_stack.pop_back();
            } catch (GRBException &e) {
                std::cerr << e.getMessage() << std::endl;
                throw e;
            }
        }

        void after_unmark(VertexPair uv) override {
            try {
                relax_pair(uv);
            } catch (GRBException &e) {
                std::cerr << e.getMessage() << std::endl;
                throw e;
            }
        }

        Cost calculate_lower_bound(Cost /*k*/) override {
            try {
                if (shall_solve) {
                    size_t result = solve();
                    if (result > initial_k) {
                        return result;
                    }
                    shall_solve = false;
                }
            } catch (GRBException &e) {
                std::cerr << e.getMessage() << std::endl;
                throw e;
            }

            return 0;
        }

/*    private:
        struct forbidden_count {
            VertexPair node_pair;
            double integrality_gap;
            double graph_diff;
            size_t num_forbidden;

            forbidden_count(VertexPair pair, double integrality_gap, double graph_diff,
                            size_t num_forbidden) : node_pair(pair), integrality_gap(integrality_gap),
                                                    graph_diff(graph_diff), num_forbidden(num_forbidden) {}

            bool operator<(const forbidden_count &other) const {
                return std::tie(this->integrality_gap, this->graph_diff, this->num_forbidden) >
                       std::tie(other.integrality_gap, other.graph_diff, other.num_forbidden);
            }

            operator VertexPair() const {
                return node_pair;
            }
        };

    public:
        ProblemSet result(const Subgraph_Stats_type &subgraph_stats, size_t k, Graph const &graph,
                          Graph const &edited, Options::Tag::Selector) {
            ProblemSet problem;
            problem.found_solution = (subgraph_stats.num_subgraphs == 0);
            problem.needs_no_edit_branch = false;
            if (!problem.found_solution && k > 0 && !shall_solve) {
                assert(model->get(GRB_IntAttr_Status) == GRB_OPTIMAL);

                double max_integrality_gap = -1;
                VertexPair uuvv = {0, 0};
                size_t max_subgraphs = 0;
                double max_diff = -1;
                for (VertexPair uv : graph.vertexPairs()) {
                    auto var = variables[uv];
                    double integrality_gap = 0.5 - std::abs(0.5 - var.get(GRB_DoubleAttr_X));
                    const double diff = std::abs((graph.hasEdge(uv) ? 1.0 : 0.0) - var.get(GRB_DoubleAttr_X));
                    const size_t subgraphs = subgraph_stats.num_subgraphs_per_edge(uv);
                    if (subgraphs == 0) return false;
                    if (integrality_gap < 0.0001) integrality_gap = 0.0;
                    if (std::tie(integrality_gap, diff, subgraphs) >
                        std::tie(max_integrality_gap, max_diff, max_subgraphs)) {
                        uuvv = uv;
                        max_integrality_gap = integrality_gap;
                        max_subgraphs = subgraphs;
                        max_diff = diff;
                    }
                    return false;
                }

                std::vector<forbidden_count> best_pairs, current_pairs;

                finder->find_near(uuvv, [&](const Subgraph &fs) {
                    current_pairs.clear();

                    for (VertexPair uv : fs.vertexPairs()) {
                        double val = variables[uv].get(GRB_DoubleAttr_X);
                        if (val < 0.0001) val = 0;
                        else if (val > 0.9999)
                            val = 1.0;
                        else if (std::abs(val - 0.5) < 0.0001)
                            val = 0.5;
                        current_pairs.emplace_back(uv, 0.5 - std::abs(0.5 - val),
                                                   std::abs((graph.hasEdge(uv) ? 1.0 : 0.0) - val),
                                                   subgraph_stats.num_subgraphs_per_edge(uv);
                        return false;
                    }

                    assert(current_pairs.size() > 0);

                    std::sort(current_pairs.begin(), current_pairs.end());

                    if (best_pairs.empty()) {
                        best_pairs = current_pairs;
                    } else {
                        size_t bi = 0, ci = 0;

                        while (bi < best_pairs.size() && ci < current_pairs.size() &&
                               (!(best_pairs[bi] < current_pairs[ci]) && !(current_pairs[ci] < best_pairs[bi]))) {
                            ++bi;
                            ++ci;
                        }

                        if (ci == current_pairs.size() ||
                            (bi != best_pairs.size() && current_pairs[ci] < best_pairs[ci])) {
                            best_pairs = current_pairs;
                        }
                    }
                }

                        assert(best_pairs.size() > 0);

                for (size_t i = 0; i < best_pairs.size(); ++i) {
                    const forbidden_count &pair_count = best_pairs[i];
                    problem.vertex_pairs.emplace_back(pair_count.node_pair, (i > 0 && i + 1 < best_pairs.size() &&
                                                                             best_pairs[i - 1].num_forbidden > 1));
                }
            }

            return problem;
        }*/
    };
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_GUROBILOWERBOUND_H
