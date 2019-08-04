//
// Created by jonas on 25.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ITERATEDLOCALSEARCH_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ITERATEDLOCALSEARCH_H


#include "../interfaces/LowerBoundI.h"

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
    Graph bound_graph;
    const VertexPairMap<Cost> &costs;
    const VertexPairMap<bool> &forbidden;
public:
    explicit IteratedLocalSearch(const Graph &graph, const VertexPairMap<Cost> &costs,
                                 const VertexPairMap<bool> &forbidden, std::shared_ptr<FinderI> finder) : LowerBoundI(
            std::move(finder)), bound_graph(graph.size()), costs(costs), forbidden(forbidden) {}

    Cost result(StateI *state, Cost k) override {
        auto s = dynamic_cast<State*>(state);
        auto &bound = s->bound;
        if (s->cost <= k) {
            initialize_bound_graph(*s);
            optimize_bound(*s, k);
        }
        return s->cost;
    }

    std::unique_ptr<StateI> initialize(Cost k) override {
        auto state = std::make_unique<State>();
        finder->find(bound_graph, [&](const Subgraph &subgraph) {
            Cost min_cost = std::numeric_limits<Cost>::max();
            bool touches_bound = subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
                min_cost = std::min(min_cost, costs[uv]);
                return bound_graph.has_edge(uv);
            });
            if (!touches_bound) {
                state->bound.emplace_back(min_cost, subgraph);
                state->cost += min_cost;
                subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
                    bound_graph.set_edge(uv);
                    return false;
                });
            }
            return state->cost > k;
        });
        return state;
    }

    void before_mark_and_edit(StateI *state, VertexPair uv) override {
        auto s = dynamic_cast<State *>(state);
        auto &bound = s->bound;

        for (size_t i = 0; i < bound.size();) {
            bool has_uv = bound[i].second.for_all_vertex_pairs([&](auto xy) {
                return uv == xy;
            });

            if (has_uv) {
                s->cost -= bound[i].first;
                bound[i] = bound.back();
                bound.pop_back();
            } else {
                ++i;
            }
        }
    }

    void after_mark_and_edit(StateI *state, VertexPair uv) override {
        auto s = dynamic_cast<State *>(state);
        auto &bound = s->bound;

        initialize_bound_graph(*s);

        finder->find_near(uv, bound_graph, [&](const Subgraph &subgraph) {
            Cost min_cost = std::numeric_limits<Cost>::max();
            bool touches_bound = subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
                min_cost = std::min(min_cost, costs[uv]);
                return bound_graph.has_edge(uv);
            });
            if (!touches_bound) {
                bound.emplace_back(min_cost, subgraph);
                s->cost += min_cost;
                subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
                    bound_graph.set_edge(uv);
                    return false;
                });
            }
            return false;
        });
    }

    void before_mark(StateI *state, VertexPair uv) override { /* no op */ }

    void after_mark(StateI *state, VertexPair uv) override { after_mark_and_edit(state, uv); }

    void before_edit(StateI *state, VertexPair uv) override { /* no op */ }

    void after_edit(StateI *state, VertexPair uv) override { /* no op */ }

    void after_unmark(StateI *state, VertexPair uv) override { /* no op */ }

private:
    void initialize_bound_graph(const State &state) {
        bound_graph.clear_edges();
        for (const auto &p : state.bound) {
            p.second.for_all_unmarked_vertex_pairs(forbidden,[&](VertexPair uv) {
                bound_graph.set_edge(uv);
                return false;
            });
        }
    }

    void optimize_bound(State &state, Cost k) {
        std::mt19937 gen(state.cost + state.bound.size());

        bool improvement_found;
        bool bound_changed;
        size_t rounds_no_improvement = 0;
        do {
            improvement_found = false;
            bound_changed = false;

            for (size_t sg_i = 0; sg_i < state.bound.size(); ++sg_i) {
                bool early_exit = find_partner(state, sg_i, k, gen, improvement_found, bound_changed);
                if (early_exit) return;
            }

            rounds_no_improvement = improvement_found ? 0 : rounds_no_improvement + 1;
        } while (improvement_found || (rounds_no_improvement < 5 && bound_changed));
    }

    bool find_partner(State &state, size_t sg_i, Cost k, std::mt19937 &gen, bool &improvement_found, bool &bound_changed) {
        constexpr Cost invalid_max_cost = std::numeric_limits<Cost>::min();

        auto subgraph = state.bound[sg_i].second;

        // remove subgraph from lower bound
        Cost subgraph_cost = cost(subgraph, forbidden, costs);
        state.cost -= subgraph_cost;
        remove_from_graph(subgraph, forbidden, bound_graph);

        // TODO: Group neighbors by vertex pairs. Subgraphs sharing a vertex pair can never be inserted together.
        // candidates are subgraphs which are only adjacent to subgraph but no other subgraph in the lower bound.
        const auto candidates = get_candidates(finder.get(), subgraph, forbidden, bound_graph);

        Cost max_subgraphs_cost = subgraph_cost;
        std::vector<size_t> max_subgraphs;

        // for each candidate check if
        //   1. the candidate has a larger cost than the current maximum cost or
        //   2. a pair of two candidates can both be inserted and their cost is larger than the current maximum cost.
        for (size_t a_i = 0; a_i < candidates.size(); ++a_i) {
            const auto &subgraph_a = candidates[a_i];

            insert_into_graph(subgraph_a, forbidden, bound_graph);

            Cost max_cost_b = invalid_max_cost;
            size_t max_b = -1;
            // for each candidate pair (a, b)
            for (size_t b_i = a_i + 1; b_i < candidates.size(); ++b_i) {
                const auto &subgraph_b = candidates[b_i];
                bool inserted = try_insert_into_graph(subgraph_b, forbidden, bound_graph);

                if (inserted) {
                    Cost cost_b = cost(subgraph_b, forbidden, costs);
                    if (cost_b > max_cost_b) {
                        max_cost_b = cost_b;
                        max_b = b_i;
                    }
                    remove_from_graph(subgraph_b, forbidden, bound_graph);
                }
            }

            Cost cost_a = cost(subgraph_a, forbidden, costs);
            if (max_b == -1) {
                // no partner found
                if (cost_a > max_subgraphs_cost) {
                    max_subgraphs_cost = cost_a;
                    max_subgraphs = {a_i};
                    improvement_found = true;
                    bound_changed = true;
                } else if (cost_a == max_subgraphs_cost) {
                    // TODO: plateau search
                }
            } else {
                // partner found
                if (cost_a + max_cost_b > max_subgraphs_cost) {
                    max_subgraphs_cost = cost_a + max_cost_b;
                    max_subgraphs = {a_i, max_b};
                    improvement_found = true;
                    bound_changed = true;
                }
            }
            remove_from_graph(subgraph_a, forbidden, bound_graph);
        }

        if (!max_subgraphs.empty()) {
            // better candidates found
            const Cost cost_0 = cost(candidates[max_subgraphs[0]], forbidden, costs);
            state.bound[sg_i] = {cost_0, candidates[max_subgraphs[0]]};
            for (size_t i = 1; i < max_subgraphs.size(); ++i) {
                const Cost cost_i = cost(candidates[max_subgraphs[i]], forbidden, costs);
                state.bound.emplace_back(cost_i, candidates[max_subgraphs[i]]);
            }
            state.cost += max_subgraphs_cost;
            if (k < state.cost) return true; // lower bound is already large enough
        } else {
            // subgraph is the best
            insert_into_graph(subgraph, forbidden, bound_graph);
            state.cost += subgraph_cost; // k stays the same
        }

        return false; // normal exit
    }

    static Cost cost(const Subgraph &subgraph, const VertexPairMap<bool> &forbidden, const VertexPairMap<Cost> &costs) {
        Cost min_vertex_pair_cost = std::numeric_limits<Cost>::max();
        subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
            min_vertex_pair_cost = std::min(min_vertex_pair_cost, costs[uv]);
        });
        assert(min_vertex_pair_cost != std::numeric_limits<Cost>::max());
        return min_vertex_pair_cost;
    }

    static void insert_into_graph(const Subgraph &subgraph, const VertexPairMap<bool>& forbidden, Graph& graph) {
        subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
            assert(!graph.has_edge(uv));
            graph.set_edge(uv);
            return false;
        });
    }

    static void remove_from_graph(const Subgraph &subgraph, const VertexPairMap<bool> &forbidden, Graph& graph) {
        subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
            assert(graph.has_edge(uv));
            graph.clear_edge(uv);
            return false;
        });
    }

    static bool try_insert_into_graph(const Subgraph& subgraph, const VertexPairMap<bool>& forbidden, Graph& graph) {
        int pairs_inserted = 0;
        bool failed = subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
            if (graph.has_edge(uv)) return true;
            graph.set_edge(uv);
            pairs_inserted++;
            return false;
        });
        if (failed) {
            subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
                assert(graph.has_edge(uv));
                graph.clear_edge(uv);
                pairs_inserted--;
                return pairs_inserted == 0;
            });
        }
        return failed;
    }

    static std::vector<Subgraph> get_candidates(FinderI* finder, const Subgraph &subgraph, const VertexPairMap<bool> &forbidden, Graph& bound_graph) {
        // Find neighboring subgraphs
        // Precondition: subgraph is removed from bound_graph.
        std::vector<Subgraph> result;

        subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
            assert(!bound_graph.has_edge(uv));
            finder->find_near(uv, bound_graph, [&](Subgraph&& neighbor) {
                result.push_back(std::move(neighbor));
                return false;
            });

            // prevent subgraphs including uv to be counted twice
            bound_graph.set_edge(uv);
        });

        // reset bound_graph
        subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
            bound_graph.clear_edge(uv);
        });

        return result;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ITERATEDLOCALSEARCH_H
