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
            std::move(finder)), bound_graph(graph.n_vertices()), costs(costs), forbidden(forbidden) {}

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
        /*
        bool improvement_found;
        bool bound_changed;
        size_t rounds_no_improvement = 0;
        do {
            improvement_found = false;
            bound_changed = false;

            for (size_t sg_i = 0; sg_i < state.bound.size(); ++sg_i) {
                bool early_exit = find_partner(state, sg_i, k, improvement_found, bound_changed);
                if (early_exit) return;
            }

            rounds_no_improvement = improvement_found ? 0 : rounds_no_improvement + 1;
        } while (improvement_found || (rounds_no_improvement < 5 && bound_changed)); */
    }
    /*
    bool find_partner(State &state, size_t sg_i, Cost k, bool &improvement_found, bool &bound_changed) {
        auto fs = state.bound[sg_i].second;

        // See how many candidates we have
        // size_t num_pairs = 0, num_neighbors = 0;
        // const auto [num_pairs, num_neighbors] = count_neighbors(subgraph_stats, g, e, fs);
        const auto pairs = get_vertex_pairs(fs, forbidden);
        // const auto [num_pairs, num_neighbors, pairs] = get_pairs(subgraph_stats, g, e, fs);

        // Skip if this forbidden subgraph if it is the only one we can add.
        // if (num_pairs > 1) {
            // Remove fs from lower bound
            remove_from_bound_uses(pairs, bound_uses);
            add_to_candidate_pairs_used(pairs, candidate_pairs_used);

            // Collect candidates
            const auto candidates_per_pair = get_candidates(finder, g, e, candidate_pairs_used, pairs, bound_uses);

            // Remove fs again from lower bound
            remove_from_bound_uses(pairs, bound_uses);
            remove_from_candidate_pairs_used(pairs, candidate_pairs_used);

            const bool random_switch = prob(gen) < 0.3;
            size_t min_candidate_neighbors = num_neighbors;
            size_t num_candidates_considered = 0;

            auto min_candidate = fs;
            size_t min_pairs = num_pairs;

            bool found_partner = false;

            for (size_t pi = 0; pi < pairs.size(); ++pi) {
                for (auto cand_fs : candidates_per_pair[pi]) {
                    assert(!found_partner);

                    add_to_bound_uses(g, e, cand_fs, bound_uses);
                    const auto [cand_pairs, cand_neighbors] = count_neighbors(subgraph_stats, g, e, cand_fs);

                    ++num_candidates_considered;

                    if (cand_pairs == 1 || (min_pairs > 1 && (
                            (!random_switch && cand_neighbors < min_candidate_neighbors) ||
                            (random_switch && prob(gen) < 1.0 / num_candidates_considered)))) {
                        min_pairs = cand_pairs;
                        min_candidate = cand_fs;
                        min_candidate_neighbors = cand_neighbors;
                    }

                    // Look only for later pairs - combinations with subgraphs listed for earlier pairs have already been considered!
                    for (size_t ppi = pi + 1; ppi < pairs.size(); ++ppi) {
                        const auto partner_pair = pairs[ppi];

                        if (bound_uses.has_edge(partner_pair.first, partner_pair.second)) continue;

                        for (auto partner_fs : candidates_per_pair[ppi]) {
                            bool touches_bound = is_subgraph_touching_bound(g, e, partner_fs, bound_uses);

                            if (!touches_bound) {
                                // Success!
                                found_partner = true;
                                improvement_found = true;

                                // Directly add the partner to the lower bound, continue search to see if there is more than one partner
                                lb.push_back(partner_fs);

                                add_to_bound_uses(g, e, partner_fs, bound_uses);
                            }
                        }

                    }

                    if (found_partner) {
                        lb[fsi] = cand_fs;

                        state.lb.assert_valid(g, e);
                        break;
                    } else {
                        remove_from_bound_uses(g, e, cand_fs, bound_uses);
                    }
                }

                if (found_partner) break;
            }

            if (!found_partner) {
                if (min_candidate != fs) {
                    add_to_bound_uses(g, e, min_candidate, bound_uses);

                    lb[fsi] = min_candidate;
                    state.lb.assert_valid(g, e);

                    bound_changed = true;
                } else {
                    // Add fs back to lower bound
                    add_to_bound_uses(pairs, bound_uses);
                }
            } else if (k < state.cost) {
                return true; //exit early
            }
        //}

        return false; // normal exit
    }

    static std::vector<VertexPair> get_vertex_pairs(const Subgraph &sg, const VertexPairMap<bool> &forbidden) {
        std::vector<VertexPair> pairs;
        sg.for_all_vertex_pairs([&](VertexPair uv) {
            if (!forbidden[uv]) pairs.push_back(uv);
            return false;
        });
        return pairs;
    }

    static void remove_from_bound_graph(const std::vector<VertexPair>& pairs, Graph &bound_graph) {
        for (auto uv : pairs) {
            bound_graph.clear_edge(uv);
        }
    }

    static std::vector<Subgraph> get_candidates(FinderI* finder, const Subgraph& sg, const VertexPairMap<bool> &forbidden, Graph &bound_graph) {
        std::vector<Subgraph> candidates;
        sg.for_all_vertex_pairs([&](VertexPair uv) {
            if (!forbidden[uv]) {
                finder->find_near(uv, bound_graph, [&](const Subgraph& neighbor) {
                    candidates.push_back(neighbor);<
                    return false;
                });
            }
            return false;
        });
        return candidates;
    }*/
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ITERATEDLOCALSEARCH_H
