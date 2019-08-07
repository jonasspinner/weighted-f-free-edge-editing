//
// Created by jonas on 25.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ITERATEDLOCALSEARCH_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ITERATEDLOCALSEARCH_H


#include "../interfaces/LowerBoundI.h"
#include "../interfaces/FinderI.h"

class IteratedLocalSearch : public LowerBoundI {
    class State : public StateI {
    public:
        std::vector<std::pair<Cost, Subgraph>> bound;
        Cost cost = 0;

        std::unique_ptr<StateI> copy() override {
            return std::make_unique<State>(*this);
        }
    };

private:
    const Graph& m_graph;
    Graph m_bound_graph;
    const VertexPairMap<Cost> &m_costs;
    const VertexPairMap<bool> &m_marked;
public:
    explicit IteratedLocalSearch(const Instance &instance,
                                 const VertexPairMap<bool> &forbidden, std::shared_ptr<FinderI> finder) : LowerBoundI(
            std::move(finder)), m_graph(instance.graph), m_bound_graph(instance.graph.size()), m_costs(instance.costs), m_marked(forbidden) {}

    /**
     *
     * @param _state
     * @param k
     * @return A lower bound on the costs required to solve the current instance.
     */
    Cost result(StateI &_state, Cost k) override {
        auto state = dynamic_cast<State &>(_state);

        state.cost = 0;
        for (size_t i = 0; i < state.bound.size(); ++i) {
            Cost min_cost = cost(state.bound[i].second, m_marked, m_costs);
            state.bound[i].first = min_cost;
            if (min_cost != invalid_cost && state.cost != invalid_cost)
                state.cost += state.bound[i].first;
            else if (min_cost == invalid_cost)
                state.cost = invalid_cost;
        }

        if (state.cost <= k) {
            initialize_bound_graph(state);
            optimize_bound(state, k);
        }
        return state.cost;
    }

    /**
     * Initializes the state by greedily constructing a maximal lower bound.
     *
     * @param k
     * @return
     */
    std::unique_ptr<StateI> initialize(Cost k) override {
        auto ptr = std::make_unique<State>();
        State& state = *ptr;

        m_bound_graph.clear_edges();

        std::vector<std::pair<Cost, Subgraph>> subgraphs;

        finder->find([&](Subgraph &&subgraph) {
            Cost min_cost = cost(subgraph, m_marked, m_costs);
            subgraphs.emplace_back(min_cost, std::move(subgraph));
            return false;
        });

        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

        for (size_t i = 0; i < subgraphs.size(); ++i) {
            bool inserted = try_insert_into_graph(subgraphs[i].second, m_marked, m_bound_graph);
            if (inserted) {
                state.cost += subgraphs[i].first;
                state.bound.push_back(std::move(subgraphs[i]));
            }
        }

        std::cout << "initial lower bound " << state.cost << "\n";

        return ptr;
    }

    /**
     * Removes subgraphs from state.bound which have the vertex pair uv.
     *
     * @param _state
     * @param uv
     */
    void before_mark_and_edit(StateI &_state, VertexPair uv) override {
        auto state = dynamic_cast<State &>(_state);
        auto &bound = state.bound;

        for (size_t i = 0; i < bound.size();) {
            const Subgraph &subgraph = bound[i].second;

            bool has_u = false; for (Vertex x : subgraph) if (x == uv.u) { has_u = true; break; }
            bool has_v = false; for (Vertex x : subgraph) if (x == uv.v) { has_v = true; break; }

            if (has_u && has_v) {
                // remove subgraph bound[i]
                state.cost -= bound[i].first;
                bound[i] = std::move(bound.back());
                bound.pop_back();
            } else {
                ++i;
            }
        }
    }

    void after_mark_and_edit(StateI &_state, VertexPair uv) override {
        auto state = dynamic_cast<State &>(_state);
        auto &bound = state.bound;

        initialize_bound_graph(state);

        std::vector<std::pair<Cost, Subgraph>> subgraphs;

        // The finder iterates over subgraphs having u and v as vertices.
        finder->find_near(uv, m_bound_graph, [&](Subgraph &&subgraph) {
            Cost min_cost = cost(subgraph, m_marked, m_costs);
            subgraphs.emplace_back(min_cost, std::move(subgraph));
            return false;
        });

        std::sort(subgraphs.begin(), subgraphs.end(), [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

        for (auto& p : subgraphs) {
            Cost min_cost = p.first;
            Subgraph& subgraph = p.second;

            // In the for loop a subgraph can be invalidated when bound_graph is modified.
            // Therefore it has to checked whether subgraph is still valid and does not share a vertex pair with bound_graph.
            bool touched = p.second.for_all_unmarked_vertex_pairs(m_marked, [&](VertexPair uv) {
                return m_bound_graph.has_edge(uv);
            });

            if (!touched) {
                insert_into_graph(subgraph, m_marked, m_bound_graph);
                bound.emplace_back(std::move(p));
                state.cost += min_cost;
            }
        }
    }

    void before_mark(StateI &state, VertexPair uv) override { /* no op */ }

    void after_mark(StateI &state, VertexPair uv) override { /* no op */ }

    void before_edit(StateI &state, VertexPair uv) override { /* no op */ }

    void after_edit(StateI &state, VertexPair uv) override { /* no op */ }

    void after_unmark(StateI &state, VertexPair uv) override { /* no op */ }

private:
    void initialize_bound_graph(const State &state) {
        m_bound_graph.clear_edges();

        for (const auto &p : state.bound) {
            p.second.for_all_unmarked_vertex_pairs(m_marked, [&](VertexPair uv) {
                m_bound_graph.set_edge(uv);
                return false;
            });
        }
    }

    void optimize_bound(State &state, Cost k) {
        std::mt19937 gen(state.cost + state.bound.size());

        auto end = [&]() {
            std::cout << "exited optimize_bound [" << state.cost << "]:\n";
            for (const auto &[c, sg] : state.bound) {
                std::cout << "\t" << sg << "\t " << c << " {";
                sg.for_all_vertex_pairs([&](VertexPair uv){
                    std::cout << " (" << uv << " e=" << m_graph.has_edge(uv) << " m=" << m_marked[uv] << " " << m_costs[uv] << ")";
                    return false;
                });
                std::cout << " }\n";
            }
        };

        bool improvement_found;
        bool bound_changed;
        size_t rounds_no_improvement = 0;
        do {
            improvement_found = false;
            bound_changed = false;

            for (size_t sg_i = 0; sg_i < state.bound.size(); ++sg_i) {
                find_2_improvement(state, sg_i, k, gen, improvement_found, bound_changed);
                if (state.cost > k) { end(); return; }
            }

            rounds_no_improvement = improvement_found ? 0 : rounds_no_improvement + 1;
        } while (improvement_found || (rounds_no_improvement < 5 && bound_changed));

        end();
    }

    /**
     * Find a 2 improvement for subgraph at index.
     *
     * @param state
     * @param index
     * @param k
     * @param gen
     * @param improvement_found
     * @param bound_changed
     * @return
     */
    void find_2_improvement(State &state, size_t index, Cost k, std::mt19937 &gen, bool &improvement_found, bool &bound_changed) {
        constexpr Cost invalid_max_cost = std::numeric_limits<Cost>::min();
        constexpr size_t invalid_index = std::numeric_limits<size_t>::max();

        const Subgraph &subgraph = state.bound[index].second;

        // remove subgraph from lower bound
        const Cost subgraph_cost = cost(subgraph, m_marked, m_costs);
        state.cost -= subgraph_cost;
        remove_from_graph(subgraph, m_marked, m_bound_graph);

        // candidates are subgraphs which are only adjacent to subgraph but no other subgraph in the lower bound.
        const auto pairs = get_pairs(subgraph, m_marked);
        assert(!subgraph.for_all_unmarked_vertex_pairs(m_marked, [&](VertexPair uv) { return m_bound_graph.has_edge(uv); }));
        auto [candidates, border] = get_candidates(*finder, pairs, m_bound_graph);

        // std::cout << "candidates:";
        std::vector<Cost> candidate_costs(candidates.size());
        for (size_t i = 0; i < candidates.size(); ++i) {
            candidate_costs[i] = cost(candidates[i], m_marked, m_costs);
        //    bool touches = candidates[i].for_all_unmarked_vertex_pairs(m_marked, [&](VertexPair uv) { return m_bound_graph.has_edge(uv); });
        //    std::cout << " (" << candidates[i] << ", " << candidate_costs[i] << ", " << touches << ")";
        }
        // std::cout << "\n";

        Cost max_subgraphs_cost = subgraph_cost;
        std::pair<size_t, size_t> max_subgraphs{invalid_index, invalid_index};

        std::vector<size_t> plateau_candidates;

        // for each candidate check if
        //   1. the candidate has a larger cost than the current maximum cost or
        //   2. a pair of two candidates can both be inserted and their cost is larger than the current maximum cost.
        for (size_t pair_i = 0; pair_i < pairs.size(); ++pair_i) {
            for (size_t a_i = border[pair_i]; a_i < border[pair_i+1]; ++a_i) {
                const Subgraph &a = candidates[a_i];

                insert_into_graph(a, m_marked, m_bound_graph);

                Cost max_cost_b = invalid_max_cost;
                size_t max_b_i = invalid_index;

                // for each candidate pair (a, b)
                for (size_t pair_j = pair_i + 1; pair_j < pairs.size(); ++pair_j) {
                    for (size_t b_i = border[pair_j]; b_i < border[pair_j + 1]; ++b_i) {
                        const Subgraph &b = candidates[b_i];

                        bool inserted = try_insert_into_graph(b, m_marked, m_bound_graph);

                        if (inserted) {
                            Cost cost_b = candidate_costs[b_i];
                            if (cost_b > max_cost_b) { max_cost_b = cost_b; max_b_i = b_i; }
                            remove_from_graph(b, m_marked, m_bound_graph);
                        }
                    }
                }

                Cost cost_a = candidate_costs[a_i];
                if (max_b_i == invalid_index) {
                    // no partner found
                    if (cost_a > max_subgraphs_cost) { max_subgraphs_cost = cost_a; max_subgraphs = {a_i, invalid_index};
                    } else if (cost_a == max_subgraphs_cost) {
                        // TODO: plateau search
                        plateau_candidates.push_back(a_i);
                    }
                } else {
                    // partner found
                    if (cost_a + max_cost_b > max_subgraphs_cost) { max_subgraphs_cost = cost_a + max_cost_b; max_subgraphs = {a_i, max_b_i}; }
                }
                remove_from_graph(a, m_marked, m_bound_graph);
            }
        }

        auto[a_i, b_i] = max_subgraphs;
        if (a_i != invalid_index) {
            improvement_found = true;
            bound_changed = true;

            // better candidates found
            std::cout << "replaced " << subgraph << " with " << candidates[a_i];
            insert_into_graph(candidates[a_i], m_marked, m_bound_graph);
            state.bound[index] = {candidate_costs[a_i], std::move(candidates[a_i])};

            if (b_i != invalid_index) {
                std::cout << " and " << candidates[b_i];
                insert_into_graph(candidates[b_i], m_marked, m_bound_graph);
                state.bound.emplace_back(candidate_costs[b_i], std::move(candidates[b_i]));
            }

            std::cout << "\n";
            state.cost += max_subgraphs_cost;
        } else {
            // subgraph is the best
            insert_into_graph(subgraph, m_marked, m_bound_graph);
            state.cost += subgraph_cost; // k stays the same
        }
    }

    /**
     * Insert unmarked vertex pairs of "subgraph" into "graph".
     *
     * @param subgraph
     * @param forbidden
     * @param graph
     */
    static void insert_into_graph(const Subgraph &subgraph, const VertexPairMap<bool> &marked, Graph &graph) {
        subgraph.for_all_unmarked_vertex_pairs(marked, [&](VertexPair uv) {
            assert(!graph.has_edge(uv));
            graph.set_edge(uv);
            return false;
        });
    }

    /**
     * Remove unmarked vertex pairs of "subgraph" from "graph".
     *
     * @param subgraph
     * @param forbidden
     * @param graph
     */
    static void remove_from_graph(const Subgraph &subgraph, const VertexPairMap<bool> &marked, Graph &graph) {
        subgraph.for_all_unmarked_vertex_pairs(marked, [&](VertexPair uv) {
            assert(graph.has_edge(uv));
            graph.clear_edge(uv);
            return false;
        });
    }

    /**
     * If possible insert all non forbidden vertex pairs of "subgraph" into "graph". If one vertex pair is already in
     * "graph", all previous insertions are inverted.
     * Returns true if all vertex pairs were successfully inserted, false otherwise.
     *
     * @param subgraph
     * @param forbidden
     * @param graph
     * @return
     */
    static bool try_insert_into_graph(const Subgraph &subgraph, const VertexPairMap<bool> &forbidden, Graph &graph) {
        // TODO: The commented code should have the same result but be more efficient. Currently the assertion is triggered.
        /* int pairs_inserted = 0;
        bool failed = subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
            if (graph.has_edge(uv)) return true;
            graph.set_edge(uv);
            pairs_inserted++;
            return false;
        });

        if (failed) {
            subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
                if (pairs_inserted == 0) return true;
                assert(graph.has_edge(uv));
                graph.clear_edge(uv);
                pairs_inserted--;
                return false;
            });
        }
        return !failed; */
        bool touches = subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) { return graph.has_edge(uv); });
        if (!touches) {
            subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) { graph.set_edge(uv); return false; });
            return true;
        } else {
            return false;
        }
    }

    /**
     * Get candidate forbidden subgraphs adjacent to a pair in "pairs".
     * Border is a separation of the candidates by the pairs.
     *
     * @param finder
     * @param pairs
     * @param bound_graph
     * @return candidates and a border array partitioning the candidates by vertex pairs.
     */
    static std::pair<std::vector<Subgraph>, std::vector<size_t>>
    get_candidates(FinderI &finder, const std::vector<VertexPair> &pairs, Graph &bound_graph) {
        // Precondition: subgraph is removed from bound_graph.
        std::vector<Subgraph> candidates;
        std::vector<size_t> border(pairs.size() + 1);

        for (size_t i = 0; i < pairs.size(); ++i) {
            VertexPair uv = pairs[i];
            assert(!bound_graph.has_edge(uv));

            finder.find_near(uv, bound_graph, [&](Subgraph &&neighbor) {
                // TODO: The condition should always be true.
                if (!neighbor.for_all_vertex_pairs([&](VertexPair uv) { return bound_graph.has_edge(uv); }))
                    candidates.push_back(std::move(neighbor));
                return false;
            });
            border[i+1] = candidates.size();

            // prevent subgraphs including uv to be counted twice
            bound_graph.set_edge(uv);
        }

        // reset bound_graph
        for (VertexPair uv : pairs) {
            bound_graph.clear_edge(uv);
        }

        return {std::move(candidates), std::move(border)};
    }

    /**
     * Return editable vertex pairs of "subgraph".
     *
     * @param subgraph
     * @param forbidden
     * @return
     */
    static std::vector<VertexPair> get_pairs(const Subgraph& subgraph, const VertexPairMap<bool> &forbidden) {
        std::vector<VertexPair> pairs;
        subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
            pairs.push_back(uv);
            return false;
        });
        return pairs;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ITERATEDLOCALSEARCH_H
