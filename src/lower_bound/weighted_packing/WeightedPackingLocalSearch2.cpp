#include "WeightedPackingLocalSearch2.h"

#include <chrono>


namespace lower_bound {

    /**
     * Returns a lower bound on the editing cost. The state keeps a maximal lower bound. Local search tries to improve
     * bound by performing (1, 1), (1, 2) or (\omega, 1) swaps on the forbidden subgraphs to improve the lower bound.
     *
     * @param k Remaining editing cost.
     * @return A lower bound on the costs required to solve the current instance.
     */
    Cost WeightedPackingLocalSearch2::calculate_lower_bound(Cost k) {
        auto &state = current_state();
        if (!state.solvable()) return state.cost();


        state.recalculate(*finder, m_marked, m_costs);
        initialize_bound_graph(*finder, state, m_marked, m_weighted_packing);
        assert(state_is_valid(*finder, state, m_marked, m_costs));
        assert(bound_graph_is_valid(*finder, state, m_marked, m_weighted_packing));

        if (state.cost() <= k)
            local_search(state, k);

        assert(state.cost() >= 0);

        return state.cost();
    }

    /**
     * Initializes the state by greedily constructing a maximal lower bound.
     *
     * @return
     */
    void WeightedPackingLocalSearch2::initialize(Cost /*k*/) {
        using Element = State::Element;

        m_states.clear();
        m_states.push_back(std::make_unique<State>());
        State &state = *m_states.back();

        // m_bound_graph.clearEdges();
        m_weighted_packing.clear();


        std::vector<Element> subgraphs;
        bool unsolveable = finder->find(m_graph, [&](Subgraph &&subgraph) {

            Cost min_cost = finder->calculate_min_cost(subgraph, m_marked, m_costs);
            subgraphs.push_back({min_cost, std::move(subgraph)});

            // If a subgraph exists with fully marked vertex pairs, the instance is not solvable.
            return min_cost == invalid_cost;
        });

        if (unsolveable) {
            state.set_unsolvable();
            return;
        }

        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const Element &lhs, const Element &rhs) { return lhs.cost > rhs.cost; });

        for (auto&[cost, subgraph] : subgraphs) {
            bool inserted = try_insert_into_packing(*finder, subgraph, m_marked, m_weighted_packing);

            if (inserted) {
                assert(cost != invalid_cost);
                state.insert({cost, std::move(subgraph)});
            }
        }

        assert(bound_graph_is_valid(*finder, state, m_marked, m_weighted_packing));
        assert(bound_is_maximal(*finder, m_graph, m_weighted_packing));
    }

    bool WeightedPackingLocalSearch2::bound_graph_is_valid(const FinderI &finder, State &state, const VertexPairMap<bool> &marked,
                                           const WeightedPacking &weighted_packing) {
        if (!state.solvable())
            return true;

        bool valid = true;
        for (VertexPair xy : Graph::VertexPairs(marked.size())) {
            if (marked[xy]) {
                if (weighted_packing.is_depleted(xy)) {
                    valid = false;
                    std::cerr << xy << " is marked and in bound_graph\n";
                }
            } else {
                /*
                bool has_subgraph = false;
                for (const auto &[cost, subgraph] : state.bound())
                    if (subgraph.contains(xy)) {
                        has_subgraph = true;
                        break;
                    }

                if (bound_graph.hasEdge(xy) != has_subgraph) {
                    valid = false;
                    if (bound_graph.hasEdge(xy))
                        std::cerr << xy << " is in bound_graph but has no subgraph in bound\n";
                    else
                        std::cerr << xy << " is not in bound_graph but has subgraph in bound\n";
                }
                */
            }
        }
        for (const auto &[cost, subgraph] : state.bound())
            finder.for_all_conversionless_edits(subgraph, [&](auto xy) {
                if (!marked[xy])
                    if (!bound_graph.hasEdge(xy)) {
                        valid = false;
                        std::cerr << xy << " is unmarked vertex pair in bound but not in bound_graph\n";
                    }
                return false;
            });
        return valid;
    }

    bool
    WeightedPackingLocalSearch2::state_is_valid(const FinderI &finder, State &state, const VertexPairMap<bool> &marked, const VertexPairMap<Cost> &costs) {
        if (!state.solvable()) return true;

        bool valid = true;

        Cost sum = 0;
        for (const auto &e : state.bound()) {
            if (e.cost != finder.calculate_min_cost(e.subgraph, marked, costs)) {
                valid = false;
            }
            sum += e.cost;
        }
        if (state.cost() != sum)
            valid = false;

        return valid;
    }

    bool WeightedPackingLocalSearch2::bound_is_maximal(FinderI &finder, const Graph &graph, const WeightedPacking &weighted_packing) {
        std::vector<Subgraph> subgraphs;
        finder.find(graph, bound_graph, [&](Subgraph &&subgraph) {
            subgraphs.push_back(std::move(subgraph));
            return false;
        });
#ifndef NDEBUG
        if (!subgraphs.empty()) {
            std::cerr << "bound is not maximal";
            for (const auto &subgraph : subgraphs)
                std::cerr << " " << subgraph;
            std::cerr << "\n";
        }
#endif
        return subgraphs.empty();
    }

    void WeightedPackingLocalSearch2::remove_near_subgraphs_from_bound(State &state, VertexPair uv) {
        assert(state.solvable());
        for (size_t i = 0; i < state.bound().size();) {
            const auto &[cost, subgraph] = state.bound(i);
            if (subgraph.contains(uv)) {
                state.remove(i);
            } else {
                ++i;
            }
        }

#ifndef NDEBUG
        for (const auto &[cost, subgraph] : state.bound())
            assert(!subgraph.contains(uv));
#endif
    }


    /**
     * Removes the subgraph from the bound containing u and v, if it exists.
     *
     * Finds forbidden subgraphs which share a vertex pair with the removed subgraph, but do not contain uv.
     *
     * @param state
     * @param uv
     */
    void WeightedPackingLocalSearch2::update_near_subgraphs(State &state, VertexPair uv, FinderI &finder,
                                            const VertexPairMap<bool> &marked, const VertexPairMap<Cost> &costs,
                                            const Graph &graph, WeightedPacking &weighted_packing) {
        assert(state.solvable());

        Subgraph removed_subgraph{};
        for (size_t i = 0; i < state.bound().size(); ++i) {
            const auto &[cost, subgraph] = state.bound(i);
            if (subgraph.contains(uv)) {
                removed_subgraph = subgraph;
                state.remove(i);
                break;
            }
        }

        if (removed_subgraph.size() == 0)
            return;

        std::vector<std::pair<Cost, Subgraph>> subgraphs;


        finder.for_all_conversionless_edits(removed_subgraph, [&](auto xy) {
            bound_graph.clearEdge(xy);
            return false;
        });
        bound_graph.setEdge(uv);

        bool unsolvable = finder.for_all_conversionless_edits(removed_subgraph, [&](auto xy) {
            if (xy == uv) return false;
            assert(!bound_graph.hasEdge(xy));

            bool unsolvable = finder.find_near_with_duplicates(xy, graph, bound_graph, [&](Subgraph &&subgraph) {  // TODO: Fix shadowing
                Cost min_cost = finder.calculate_min_cost(subgraph, marked, costs);
                subgraphs.emplace_back(min_cost, std::move(subgraph));
                return min_cost == invalid_cost;
            });

            if (unsolvable) {
                state.set_unsolvable();
                return true;
            }

            weighted_packing.setEdge(xy);
            return false;
        });

        if (unsolvable)
            return;

        finder.for_all_conversionless_edits(removed_subgraph, [&](auto xy) {
            weighted_packing.clearEdge(xy);
            return false;
        });

#ifndef NDEBUG  // TODO: Adapt to conversionless edit optimization.
//        for (const auto &[cost, subgraph] : subgraphs)
//            assert(!subgraph.contains(uv));
#endif

        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

        insert_subgraphs_into_packing(finder, std::move(subgraphs), marked, state, weighted_packing);

#ifndef NDEBUG  // TODO: Adapt to conversionless edit optimization.
//        for (const auto &[cost, subgraph] : state.bound())
//            assert(!subgraph.contains(uv)); // Throws in test
#endif
    }

    /**
     * Finds forbidden subgraph containing u and v. Tries to insert them into the bound in descending order by minimum
     * editing cost.
     *
     * If a fully marked subgraph is found, the state is marked as unsolvable.
     */
    void WeightedPackingLocalSearch2::insert_near_subgraphs_into_bound(State &state, VertexPair uv, FinderI &finder,
                                                       const VertexPairMap<bool> &marked,
                                                       const VertexPairMap<Cost> &costs, const Graph &graph,
                                                       WeightedPacking &weighted_packing) {
        assert(state.solvable());

        std::vector<std::pair<Cost, Subgraph>> subgraphs;

        // The finder iterates over subgraphs having u and v as vertices.
        bool unsolvable = finder.find_near_with_duplicates(uv, graph, weighted_packing.m_depleted_graph, [&](Subgraph &&subgraph) {
            Cost min_cost = finder.calculate_min_cost(subgraph, marked, costs);  // Only consider conversionless edits.
            subgraphs.emplace_back(min_cost, std::move(subgraph));
            return min_cost == invalid_cost;
        });

        if (unsolvable) {
            state.set_unsolvable();
            return;
        }

        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

        for (auto &&[cost, subgraph] : subgraphs) {
            // In the for loop a subgraph can be invalidated when bound_graph is modified.
            // Therefore it has to checked whether subgraph is still valid and does not share a vertex pair with bound_graph.
            bool inserted = try_insert_into_packing(finder, subgraph, marked, weighted_packing);
            if (inserted) {
                state.insert({cost, std::move(subgraph)});
            }
        }

#ifndef NDEBUG
        // local
        bool has_nearby_subgraph = finder.find_near_with_duplicates(uv, graph, weighted_packing.m_depleted_graph, [](Subgraph &&) { return true; });
        assert(!has_nearby_subgraph);
#endif
    }

    /**
     * Initializes the bound_graph. Every unmarked vertex pair of a subgraph in the bound becomes an edge in the graph.
     *
     * @param state
     * @param marked
     * @param bound_graph
     */
    void WeightedPackingLocalSearch2::initialize_bound_graph(const FinderI &finder, const State &state,
                                             const VertexPairMap<bool> &marked, WeightedPacking &weighted_packing) {
        weighted_packing.clear();
        for (const auto &[cost, subgraph] : state.bound())
            weighted_packing.subtract_from_potential(subgraph, cost);
    }

    /**
     * When a subgraph containing u and v exists, remove it from the lower bound.
     *
     * if \exists S \in B: u, v \in S then
     *      remove S from B
     *
     * @param uv
     */
    void WeightedPackingLocalSearch2::before_mark(VertexPair uv) {
        auto &state = current_state();

        if (!state.solvable()) return;

        // Either only remove the subgraph at uv (I) or remove it and fill the "hole" with subgraphs not containing uv (II).
        // remove_near_subgraphs_from_bound(state, uv); // I
        update_near_subgraphs(state, uv, *finder, m_marked, m_costs, m_graph, m_weighted_packing); // II

#ifndef NDEBUG  // TODO: Adapt to conversionless edit optimization.
//        for (const auto &[cost, subgraph] : state.bound())
//            assert(!subgraph.contains(uv));
#endif
    }

    /**
     * Update the state one recursion level higher up. Insert forbidden subgraphs which contain u and v.
     *
     * @param uv
     */
    void WeightedPackingLocalSearch2::before_edit(VertexPair uv) {
        auto &state = parent_state();

        if (!state.solvable()) return;

        initialize_bound_graph(*finder, state, m_marked, m_weighted_packing);
        insert_near_subgraphs_into_bound(state, uv, *finder, m_marked, m_costs, m_graph, m_weighted_packing);
    }

    /**
     * Update the current state. Insert the new created forbidden subgraphs. They contain u and v.
     * @param uv
     */
    void WeightedPackingLocalSearch2::after_edit(VertexPair uv) {
        auto &state = current_state();

        if (!state.solvable()) return;

        initialize_bound_graph(*finder, state, m_marked, m_weighted_packing);
        insert_near_subgraphs_into_bound(state, uv, *finder, m_marked, m_costs, m_graph, m_weighted_packing);
    }

    /**
     * Tries to optimise the bound. Stops if the lower bound is larger than k or the bound does not improve in five consecutive rounds.
     *
     * @param state
     * @param k
     */
    void WeightedPackingLocalSearch2::local_search(State &state, Cost k) {
        state.shuffle(m_gen);

        bool improvement_found, bound_changed;
        size_t rounds_no_improvement = 0;

        do {
            // a single round consists of a loop over all subgraphs in the bound and trying to find an improvement depending on the current mode.
            improvement_found = bound_changed = false;
            for (size_t index = 0; state.cost() <= k && index < state.bound().size(); ++index)
                improvement_found |= find_two_improvement(state, index, bound_changed);

            rounds_no_improvement = improvement_found ? 0 : rounds_no_improvement + 1;
            // improvement_found || (rounds_no_improvement < 5 && bound_changed)
        } while (state.cost() <= k &&
                 (improvement_found || (rounds_no_improvement < m_max_rounds_no_improvement && bound_changed)));
        if (verbosity) std::cout << "[" << state.cost() << "] end local search\n";


#ifndef NDEBUG
        if (state.cost() <= k) {
            assert(bound_is_maximal(*finder, m_graph, m_weighted_packing));
        }
#endif
    }

    /**
     * Try to find a (1, 2) swap for the subgraph at the index which improves the lower bound.
     *
     * For the subgraph currently in the lower bound at the given index do the following. List other forbidden subgraphs
     * which are adjacent to the subgraph but not to an other subgraph in the lower bound as candidates. Test for each
     * pair of candidates whether they are adjacent.
     * + From all pairs of candidates which are not adjacent, choose the pair with the largest combined cost.
     * + If single candidate has a larger cost choose that.
     * + If no improvement can be found choose a candidate with the same cost as the original subgraph at random.
     * + If no candidate can be found, let the original subgraph remain.
     *
     * Complexity: O((p over 2) * find_near_with_duplicates + |candidates|^2)
     *
     * @param state
     * @param index
     * @param gen
     * @param bound_changed
     * @return
     */
    bool WeightedPackingLocalSearch2::find_two_improvement(State &state, size_t index, bool &bound_changed) {
        constexpr Cost invalid_max_cost = std::numeric_limits<Cost>::min();
        constexpr size_t invalid_index = std::numeric_limits<size_t>::max();

        bool found_improvement = false;

        const auto &[subgraph_cost, subgraph] = state.bound(index);
        assert(subgraph_cost == finder->calculate_min_cost(subgraph, m_marked, m_costs));
        assert(subgraph_cost != invalid_cost);

        // If the subgraph has no neighbors on at least one vertex pair it can be skipped.
        const auto[num_pairs_with_neighbors, num_neighbors_upper_bound] =
            count_neighbors(*finder, m_subgraph_stats, m_marked, subgraph);
        if (num_pairs_with_neighbors < 1 || num_neighbors_upper_bound < 1)
            return false; // No improvement possible.

        // Remove the subgraph from m_bound_graph. Either the subgraph or other subgraphs will be reinserted.
        remove_from_packing(*finder, subgraph, m_marked, m_weighted_packing);

        // Candidates are subgraphs which are only adjacent to subgraph but no other subgraph in the lower bound.
        const auto pairs = get_pairs(*finder, subgraph, m_marked);
        auto[candidates, border] = get_candidates(*finder, pairs, m_graph, m_weighted_packing, m_subgraph_stats);

#ifndef NDEBUG
        {
            bool touches = finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
                return !m_marked[uv] && m_bound_graph.is_depleted(uv);
            });
            assert(!touches);
        }
#endif

        std::vector<Cost> candidate_costs(candidates.size());
        for (size_t i = 0; i < candidates.size(); ++i)
            candidate_costs[i] = finder->calculate_min_cost(candidates[i], m_marked, m_costs);


        // The information about the best solution.
        Cost max_subgraphs_cost = subgraph_cost;
        std::vector<size_t> max_subgraphs;

        std::vector<size_t> plateau_candidates;

        // for each candidate check if
        //   1. the candidate has a larger cost than the current maximum cost or
        //   2. a pair of two candidates can both be inserted and their cost is larger than the current maximum cost.
        for (size_t pair_i = 0; pair_i < pairs.size(); ++pair_i) {
            for (size_t a_i = border[pair_i]; a_i < border[pair_i + 1]; ++a_i) {
                const Subgraph &a = candidates[a_i];

                insert_into_packing(*finder, a, m_marked, m_weighted_packing);

                Cost max_cost_b = invalid_max_cost;
                size_t max_b_i = invalid_index;

                // for each candidate pair (a, b)
                for (size_t pair_j = pair_i + 1; pair_j < pairs.size(); ++pair_j) {
                    if (m_weighted_packing.is_depleted(pairs[pair_j])) continue;
                    for (size_t b_i = border[pair_j]; b_i < border[pair_j + 1]; ++b_i) {
                        const Subgraph &b = candidates[b_i];

                        bool inserted = try_insert_into_packing(*finder, b, m_marked, m_weighted_packing);

                        if (inserted) {
                            Cost cost_b = candidate_costs[b_i];
                            if (cost_b > max_cost_b) {
                                max_cost_b = cost_b;
                                max_b_i = b_i;
                            }
                            remove_from_packing(*finder, b, m_marked, m_weighted_packing);
                        }
                    }
                }

                Cost cost_a = candidate_costs[a_i];
                if (max_b_i == invalid_index) {
                    // no partner found
                    if (cost_a > max_subgraphs_cost) {
                        max_subgraphs_cost = cost_a;
                        max_subgraphs = {a_i};
                    } else if (cost_a == max_subgraphs_cost) {
                        // TODO: plateau search
                        plateau_candidates.push_back(a_i);
                    }
                } else {
                    // partner found
                    if (cost_a + max_cost_b > max_subgraphs_cost) {
                        max_subgraphs_cost = cost_a + max_cost_b;
                        max_subgraphs = {a_i, max_b_i};
                    }
                }
                remove_from_packing(*finder, a, m_marked, m_weighted_packing);
            }
        }

        if (!max_subgraphs.empty()) {
            size_t a_i = max_subgraphs[0];
            found_improvement = true;
            bound_changed = true;

            if (verbosity > 1)
                std::cout << "found (1, " << max_subgraphs.size() << ") swap " << std::setw(4)
                          << max_subgraphs_cost - subgraph_cost << ", " << subgraph_cost << " => " << max_subgraphs_cost
                          << ", " << subgraph << " => ...\n";

            // better candidates found
            insert_into_packing(*finder, candidates[a_i], m_marked, m_weighted_packing);
            state.replace(index, {candidate_costs[a_i], std::move(candidates[a_i])});

            if (max_subgraphs.size() == 2) {
                size_t b_i = max_subgraphs[1];

                insert_into_packing(*finder, candidates[b_i], m_marked, m_weighted_packing);
                state.insert({candidate_costs[b_i], std::move(candidates[b_i])});
            }

        } else {
            // subgraph is the best
            if (!plateau_candidates.empty()) {
                bound_changed = true;

                std::uniform_real_distribution<float> uniform(0, 1);

                if (uniform(m_gen) < m_alpha) {
                    // Every \alpha percent.
                    // Do not choose between all candidates, but only between the candidates with the lowest upper bound on
                    // subgraphs covered.

                    std::vector<size_t> new_plateau_candidates;
                    size_t min_num_subgraphs_covered = std::numeric_limits<size_t>::max();

                    for (auto i : plateau_candidates) {
                        auto[num_pairs, num_covered_ub] =
                            count_neighbors(*finder, m_subgraph_stats, m_marked, candidates[i]);
                        if (num_covered_ub < min_num_subgraphs_covered) {
                            new_plateau_candidates = {i};
                            min_num_subgraphs_covered = num_covered_ub;
                        } else if (num_covered_ub == min_num_subgraphs_covered) {
                            new_plateau_candidates.push_back(i);
                        }
                    }
                    plateau_candidates = new_plateau_candidates;
                }

                std::uniform_int_distribution<size_t> sample(0, plateau_candidates.size() - 1);
                auto a_i = plateau_candidates[sample(m_gen)];

                if (verbosity > 1)
                    std::cout << "made (1, 1) swap for plateau search " << std::setw(4) << 0 << ", " << subgraph
                              << " => " << candidates[a_i] << "\n";

                insert_into_packing(*finder, candidates[a_i], m_marked, m_weighted_packing);
                state.replace(index, {candidate_costs[a_i], std::move(candidates[a_i])});
            } else {
                insert_into_packing(*finder, subgraph, m_marked, m_weighted_packing);
            }
        }
        if (state.solvable())
            assert(bound_graph_is_valid(*finder, state, m_marked, m_weighted_packing));

        return found_improvement;
    }

    /**
     * Insert unmarked vertex pairs of "subgraph" into "graph".
     *
     * @param subgraph
     * @param forbidden
     * @param graph
     */
    void WeightedPackingLocalSearch2::insert_into_packing(const FinderI &finder, const Subgraph &subgraph,
                                        const VertexPairMap<bool> &marked, WeightedPacking &weighted_packing) {
        Cost cost = weighted_packing.calculate_min_cost(subgraph);
        weighted_packing.subtract_from_potential(subgraph, cost);
    }

    /**
     * Remove unmarked vertex pairs of "subgraph" from "graph".
     *
     * @param subgraph
     * @param forbidden
     * @param graph
     */
    void WeightedPackingLocalSearch2::remove_from_packing(const FinderI &finder, const Subgraph &subgraph,
                                        const VertexPairMap<bool> &marked, WeightedPacking &weighted_packing) {
        finder.for_all_conversionless_edits(subgraph, [&](auto uv) {
            if (!marked[uv]) {
                assert(graph.hasEdge(uv));
                graph.clearEdge(uv);
            }
            return false;
        });
    }

    /**
     * If possible insert all non forbidden vertex pairs of subgraph into graph. Returns true if all vertex pairs were
     * successfully inserted, false otherwise.
     *
     * @param subgraph
     * @param marked
     * @param graph
     * @return
     */
    bool WeightedPackingLocalSearch2::try_insert_into_packing(const FinderI &finder, const Subgraph &subgraph,
                                            const VertexPairMap<bool> &marked, WeightedPacking &weighted_packing) {

        Cost cost = weighted_packing.calculate_min_cost(subgraph);

        if (cost != invalid_cost) {
            weighted_packing.subtract_from_potential(subgraph, cost);
            return true;
        }
        return false;
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
    std::pair<std::vector<Subgraph>, std::vector<size_t>>
    WeightedPackingLocalSearch2::get_candidates(FinderI &finder, const std::vector<VertexPair> &pairs, const Graph &graph,
                                WeightedPacking &weighted_packing, const SubgraphStats &subgraph_stats) {

        // Precondition: subgraph is removed from bound_graph.
#ifndef NDEBUG
        for (auto uv : pairs) {
            assert(!bound_graph.hasEdge(uv));
        }
#endif

        std::vector<Subgraph> candidates;
        std::vector<size_t> border(pairs.size() + 1);

        for (size_t i = 0; i < pairs.size(); ++i) {
            VertexPair uv = pairs[i];
            assert(!bound_graph.hasEdge(uv));

            if (subgraph_stats.subgraphCount(uv) > 1) {
                finder.find_near_with_duplicates(uv, graph, bound_graph, [&](Subgraph &&neighbor) {
#ifndef NDEBUG
                    finder.for_all_conversionless_edits(neighbor, [&](auto xy) {
                        assert(!bound_graph.hasEdge(xy));
                        return false;
                    });
#endif
                    candidates.push_back(std::move(neighbor));
                    return false;
                });
            }
            border[i + 1] = candidates.size();

            // prevent subgraphs including uv to be counted twice
            bound_graph.setEdge(uv);
        }

        // reset bound_graph
        for (VertexPair uv : pairs) {
            assert(bound_graph.hasEdge(uv));
            bound_graph.clearEdge(uv);
        }

        return {std::move(candidates), std::move(border)};
    }

    /**
     * Return editable vertex pairs of subgraph.
     *
     * @param subgraph
     * @param marked
     * @return
     */
    std::vector<VertexPair> WeightedPackingLocalSearch2::get_pairs(const FinderI &finder, const Subgraph &subgraph,
                                                   const VertexPairMap<bool> &marked) {
        std::vector<VertexPair> pairs;
        finder.for_all_conversionless_edits(subgraph, [&](auto uv) {
            if (!marked[uv])
                pairs.push_back(uv);
            return false;
        });
        return pairs;
    }

    /**
     * Return the number of pairs with neighbors and a upper bound on the number of subgraphs sharing an unmarked pair
     * of vertices with the given subgraph.
     *
     * @param subgraph_stats
     * @param marked
     * @param subgraph
     * @return
     */
    std::tuple<size_t, size_t>
    WeightedPackingLocalSearch2::count_neighbors(const FinderI &finder, const SubgraphStats &subgraph_stats,
                                 const VertexPairMap<bool> &marked, const Subgraph &subgraph) {
        size_t num_pairs = 0, num_neighbors_ub = 0;
        finder.for_all_conversionless_edits(subgraph, [&](auto uv) {
            if (!marked[uv]) {
                size_t nn = subgraph_stats.subgraphCount(uv) - 1;
                num_neighbors_ub += nn;
                if (nn > 0) ++num_pairs;
            }
            return false;
        });
        return {num_pairs, num_neighbors_ub};
    }

}
