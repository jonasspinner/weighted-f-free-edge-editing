#include "LocalSearch.h"

#include <chrono>


namespace lower_bound {

    /**
     * Returns a lower bound on the editing cost. The state keeps a maximal lower bound. Local search tries to improve
     * bound by performing (1, 1), (1, 2) or (\omega, 1) swaps on the forbidden subgraphs to improve the lower bound.
     *
     * @param k Remaining editing cost.
     * @return A lower bound on the costs required to solve the current instance.
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    Cost LocalSearch<SetOfForbiddenSubgraphs>::calculate_lower_bound(Cost k) {
        auto &state = current_state();
        if (!state.solvable()) return state.cost();


        state.recalculate(m_edit_state->marked_map(), m_edit_state->cost_map());
        initialize_bound_graph(state, m_edit_state->marked_map(), m_bound_graph);
        assert(state_is_valid(state, m_edit_state->marked_map(), m_edit_state->cost_map()));
        assert(bound_graph_is_valid(state, m_edit_state->marked_map(), m_bound_graph, m_edit_state->cost_map()));

        if (state.cost() <= k)
            local_search(state, k);

        assert(state.cost() >= 0);

        return state.cost();
    }

    template<Options::FSG SetOfForbiddenSubgraphs>
    Cost LocalSearch<SetOfForbiddenSubgraphs>::calculate_lower_bound_no_edit_branch() {
        auto &state = current_state();
        return state.cost();
    }

    /**
     * Initializes the state by greedily constructing a maximal lower bound.
     *
     * @return
     */

    template<Options::FSG SetOfForbiddenSubgraphs>
    void LocalSearch<SetOfForbiddenSubgraphs>::initialize(Cost /*k*/) {

        m_states.clear();
        m_states.push_back(std::make_unique<State>());
        State &state = *m_states.back();

        m_bound_graph.clear_edges();

        for (auto uv : m_edit_state->cost_map().keys()) {
            if (m_edit_state->cost(uv) == 0) {
                m_bound_graph.set_edge(uv);
            }
        }

        std::vector<std::pair<Cost, Subgraph>> subgraphs;
        auto exit_state = m_finder.find(m_edit_state->graph(), m_bound_graph, [&](const auto &subgraph) {
            Cost min_cost = subgraph.calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map());
            subgraphs.emplace_back(min_cost, subgraph);
            return subgraph_iterators::break_if(min_cost == invalid_cost);
        });

        if (exit_state == subgraph_iterators::IterationExit::Break) {
            state.set_unsolvable();
            return;
        }

        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

        for (auto&[cost, subgraph] : subgraphs) {
            bool inserted = try_insert_into_graph(subgraph, m_edit_state->marked_map(), m_bound_graph);

            if (inserted) {
                assert(cost == subgraph.calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map()));
                assert(cost != invalid_cost);
                state.insert({cost, subgraph});
            }
        }

        assert(bound_graph_is_valid(state, m_edit_state->marked_map(), m_bound_graph, m_edit_state->cost_map()));
        assert(bound_is_maximal(m_finder, m_edit_state->graph(), m_bound_graph));
    }

    template<Options::FSG SetOfForbiddenSubgraphs>
    bool LocalSearch<SetOfForbiddenSubgraphs>::bound_graph_is_valid(State &state, const VertexPairMap<bool> &marked,
            const Graph &bound_graph, const VertexPairMap<Cost> &costs) {
        if (!state.solvable())
            return true;

        bool valid = true;
        for (VertexPair xy : Graph::VertexPairs(marked.size())) {
            if (marked[xy]) {
                if (bound_graph.has_edge(xy)) {
                    valid = false;
                    std::cerr << xy << " is marked and in bound_graph\n";
                }
            } else {
                bool has_subgraph = false;
                for (const auto &[cost, subgraph] : state.bound()) {
                    auto edits = subgraph.non_converting_edits();
                    auto contains_pair = std::any_of(edits.begin(), edits.end(), [&](auto e) {
                        return e == xy;
                    });
                    if (contains_pair) {
                        has_subgraph = true;
                        break;
                    }
                }

                if (bound_graph.has_edge(xy) != has_subgraph && costs[xy] > 0) {
                    valid = false;
                    if (bound_graph.has_edge(xy))
                        std::cerr << xy << " is in bound_graph but has no subgraph in bound\n";
                    else
                        std::cerr << xy << " is not in bound_graph but has subgraph in bound\n";
                }
            }
        }
        for (const auto &[cost, subgraph] : state.bound())
            for (auto xy : subgraph.non_converting_edits()) {
                if (!marked[xy])
                    if (!bound_graph.has_edge(xy)) {
                        valid = false;
                        std::cerr << xy << " is unmarked vertex pair in bound but not in bound_graph\n";
                    }
            }
        return valid;
    }

    template<Options::FSG SetOfForbiddenSubgraphs>
    bool LocalSearch<SetOfForbiddenSubgraphs>::state_is_valid(
            State &state, const VertexPairMap<bool> &marked, const VertexPairMap<Cost> &costs) {
        if (!state.solvable()) return true;

        bool valid = true;

        Cost sum = 0;
        for (const auto &e : state.bound()) {
            if (e.cost != e.subgraph.calculate_min_cost(costs, marked)) {
                valid = false;
            }
            sum += e.cost;
        }
        if (state.cost() != sum)
            valid = false;

        return valid;
    }

    template<Options::FSG SetOfForbiddenSubgraphs>
    bool LocalSearch<SetOfForbiddenSubgraphs>::bound_is_maximal(Finder &finder, const Graph &graph,
                                                                const Graph &bound_graph) {
        std::vector<Subgraph> subgraphs;
        finder.find(graph, bound_graph, [&](const Subgraph &subgraph) {
            subgraphs.push_back(std::move(subgraph));
            return subgraph_iterators::IterationControl::Continue;
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


    /**
     * Removes the subgraph from the bound containing u and v, if it exists.
     *
     * Finds forbidden subgraphs which share a vertex pair with the removed subgraph, but do not contain uv.
     *
     * @param state
     * @param uv
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    void LocalSearch<SetOfForbiddenSubgraphs>::update_near_subgraphs(State &state, VertexPair uv, Finder &finder,
                                                                     const VertexPairMap<bool> &marked,
                                                                     const VertexPairMap<Cost> &costs,
                                                                     const Graph &graph, Graph &bound_graph) {
        assert(state.solvable());

        std::optional<Subgraph> removed_subgraph_opt;
        for (size_t i = 0; i < state.bound().size(); ++i) {
            const auto &[cost, subgraph] = state.bound(i);
            auto edits = subgraph.non_converting_edits();
            auto contains_pair = std::any_of(edits.begin(), edits.end(), [&](auto e) {
                return e == uv;
            });
            if (contains_pair) {
                removed_subgraph_opt = subgraph;
                state.remove(i);
                break;
            }
        }

        if (!removed_subgraph_opt.has_value())
            return;

        auto removed_subgraph = *removed_subgraph_opt;

        std::vector<std::pair<Cost, Subgraph>> subgraphs;


        for (auto xy : removed_subgraph.non_converting_edits()) {
            bound_graph.reset_edge(xy);
        }
        bound_graph.set_edge(uv);

        bool unsolvable = false;
        for (auto xy : removed_subgraph.non_converting_edits()) {
            if (xy == uv)
                continue;
            assert(!bound_graph.has_edge(xy));

            auto exit_state = finder.find_near(xy, graph, bound_graph, [&](const Subgraph &subgraph) {
                Cost min_cost = subgraph.calculate_min_cost(costs, marked);
                subgraphs.emplace_back(min_cost, subgraph);
                return subgraph_iterators::break_if(min_cost == invalid_cost);
            });

            unsolvable = exit_state == subgraph_iterators::IterationExit::Break;

            if (unsolvable) {
                state.set_unsolvable();
                break;
            }

            bound_graph.set_edge(xy);
        }

        if (unsolvable)
            return;

        for (auto xy  : removed_subgraph.non_converting_edits()) {
            bound_graph.reset_edge(xy);
        }

#ifndef NDEBUG
        for (const auto &[cost, subgraph] : subgraphs) {
            auto edits = subgraph.non_converting_edits();
            auto contains_pair = std::any_of(edits.begin(), edits.end(), [&](auto e) {
                return e == uv;
            });
            assert(!contains_pair);
        }
#endif

        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

        insert_subgraphs_into_bound(std::move(subgraphs), marked, state, bound_graph);

#ifndef NDEBUG
        for (const auto &[cost, subgraph] : state.bound()) {
            auto edits = subgraph.non_converting_edits();
            auto contains_pair = std::any_of(edits.begin(), edits.end(), [&](auto e) {
                return e == uv;
            });
            assert(!contains_pair);
        }
#endif
    }

    /**
     * Finds forbidden subgraph containing u and v. Tries to insert them into the bound in descending order by minimum
     * editing cost.
     *
     * If a fully marked subgraph is found, the state is marked as unsolvable.
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    void
    LocalSearch<SetOfForbiddenSubgraphs>::insert_near_subgraphs_into_bound(State &state, VertexPair uv, Finder &finder,
                                                                           const VertexPairMap<bool> &marked,
                                                                           const VertexPairMap<Cost> &costs,
                                                                           const Graph &graph,
                                                                           Graph &bound_graph) {
        assert(state.solvable());

        std::vector<std::pair<Cost, Subgraph>> subgraphs;

        // The finder iterates over subgraphs having u and v as vertices.
        auto exit_state = finder.find_near(uv, graph, bound_graph, [&](const Subgraph &subgraph) {
            Cost min_cost = subgraph.calculate_min_cost(costs, marked);  // Only consider conversionless edits.
            subgraphs.emplace_back(min_cost, subgraph);
            return subgraph_iterators::break_if(min_cost == invalid_cost);
        });

        if (exit_state == subgraph_iterators::IterationExit::Break) {
            state.set_unsolvable();
            return;
        }

        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

        for (auto &&[cost, subgraph] : subgraphs) {
            // In the for loop a subgraph can be invalidated when bound_graph is modified.
            // Therefore it has to checked whether subgraph is still valid and does not share a vertex pair with bound_graph.
            bool inserted = try_insert_into_graph(subgraph, marked, bound_graph);
            if (inserted) {
                assert(cost == subgraph.calculate_min_cost(costs, marked));
                state.insert({cost, subgraph});
            }
        }

#ifndef NDEBUG
        // local
        bool has_nearby_subgraph =
                subgraph_iterators::IterationExit::Break == finder.find_near(uv, graph, bound_graph,
                    [](const Subgraph &) { return subgraph_iterators::IterationControl::Break; });
        assert(!has_nearby_subgraph);
#endif
    }

    /**
     * Tries to insert the subgraphs into the lower bound by the given order.
     *
     * @param subgraphs
     * @param marked
     * @param state
     * @param bound_graph
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    void LocalSearch<SetOfForbiddenSubgraphs>::insert_subgraphs_into_bound(
            std::vector<std::pair<Cost, Subgraph>> &&subgraphs,
            const VertexPairMap<bool> &marked, State &state, Graph &bound_graph) {
        for (auto &&[cost, subgraph] : subgraphs) {
            // In the for loop a subgraph can be invalidated when bound_graph is modified.
            // Therefore it has to checked whether subgraph is still valid and does not share a vertex pair with bound_graph.
            bool inserted = try_insert_into_graph(subgraph, marked, bound_graph);
            if (inserted) {
                state.insert({cost, subgraph});
            }
        }
    }

    /**
     * Initializes the bound_graph. Every unmarked vertex pair of a subgraph in the bound becomes an edge in the graph.
     *
     * @param state
     * @param marked
     * @param bound_graph
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    void LocalSearch<SetOfForbiddenSubgraphs>::initialize_bound_graph(const State &state,
                                                                      const VertexPairMap<bool> &marked,
                                                                      Graph &bound_graph) {
        bound_graph.clear_edges();
        for (const auto &[cost, subgraph] : state.bound())
            for (auto uv : subgraph.non_converting_edits()) {
                if (!marked[uv])
                    bound_graph.set_edge(uv);
            }
    }

    /**
     * When a subgraph containing u and v exists, remove it from the lower bound.
     *
     * if \exists S \in B: u, v \in S then
     *      remove S from B
     *
     * @param uv
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    void LocalSearch<SetOfForbiddenSubgraphs>::before_mark(VertexPair uv) {
        auto &state = current_state();

        if (!state.solvable()) return;

        // Either only remove the subgraph at uv (I) or remove it and fill the "hole" with subgraphs not containing uv (II).
        // remove_near_subgraphs_from_bound(state, uv); // I
        update_near_subgraphs(state, uv, m_finder, m_edit_state->marked_map(), m_edit_state->cost_map(), m_edit_state->graph(), m_bound_graph); // II

#ifndef NDEBUG
        for (const auto &[cost, subgraph] : state.bound()) {
            auto edits = subgraph.non_converting_edits();
            auto contains_pair = std::any_of(edits.begin(), edits.end(), [&](auto e) {
                return e == uv;
            });
            assert(!contains_pair);
        }
#endif
    }

    /**
     * Update the state one recursion level higher up. Insert forbidden subgraphs which contain u and v.
     *
     * @param uv
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    void LocalSearch<SetOfForbiddenSubgraphs>::before_edit(VertexPair uv) {
        auto &state = parent_state();

        if (!state.solvable()) return;

        initialize_bound_graph(state, m_edit_state->marked_map(), m_bound_graph);
        insert_near_subgraphs_into_bound(state, uv, m_finder, m_edit_state->marked_map(), m_edit_state->cost_map(), m_edit_state->graph(), m_bound_graph);
    }

    /**
     * Update the current state. Insert the new created forbidden subgraphs. They contain u and v.
     * @param uv
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    void LocalSearch<SetOfForbiddenSubgraphs>::after_edit(VertexPair uv) {
        auto &state = current_state();

        if (!state.solvable()) return;

        initialize_bound_graph(state, m_edit_state->marked_map(), m_bound_graph);
        insert_near_subgraphs_into_bound(state, uv, m_finder, m_edit_state->marked_map(), m_edit_state->cost_map(), m_edit_state->graph(), m_bound_graph);
    }

    /**
     * Tries to optimise the bound. Stops if the lower bound is larger than k or the bound does not improve in five consecutive rounds.
     *
     * @param state
     * @param k
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    void LocalSearch<SetOfForbiddenSubgraphs>::local_search(State &state, Cost k) {
        // TODO: Remove ugly debug output
// #define stats
#ifdef stats
        Cost cost_before = state.cost();
        size_t num_rounds_one = 0; Cost improvement_one = 0;
        size_t num_rounds_two = 0; Cost improvement_two = 0;
        size_t num_rounds_omega = 0; Cost improvement_omega = 0;

        long time_one = 0;
        long time_two = 0;
        long time_omega = 0;

        size_t num_marked = 0;
        for (VertexPair uv : m_bound_graph.vertexPairs())
            if (m_marked[uv])
                ++num_marked;
        std::cerr << "depth " << std::setw(2) << m_states.size() << " ";
        std::cerr << "num_marked " << std::setw(3) << num_marked << " ";
#endif

        state.shuffle(m_gen);

        bool improvement_found, bound_changed;
        size_t rounds_no_improvement = 0;

        do {
            // a single round consists of a loop over all subgraphs in the bound and trying to find an improvement depending on the current mode.
            improvement_found = bound_changed = false;


            if (use_one_improvement) {
#ifdef stats
                Cost before = state.cost();
                auto start = std::chrono::steady_clock::now();
#endif
                if (verbosity) std::cout << "[" << state.cost() << "] round mode=0\n";
                for (size_t index = 0; state.cost() <= k && index < state.bound().size(); ++index)
                    improvement_found |= find_one_improvements(state, index);
#ifdef stats
                ++num_rounds_one; improvement_one += state.cost() - before;
                time_one = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count();
#endif
            }

            if (!improvement_found && use_two_improvement) {
#ifdef stats
                Cost before = state.cost();
                auto start = std::chrono::steady_clock::now();
#endif
                if (verbosity) std::cout << "[" << state.cost() << "] round mode=1\n";
                for (size_t index = 0; state.cost() <= k && index < state.bound().size(); ++index)
                    improvement_found |= find_two_improvement(state, index, bound_changed);
#ifdef stats
                ++num_rounds_two; improvement_two += state.cost() - before;
                time_two = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count();
#endif
            }

            if (!improvement_found && use_omega_improvement) {
#ifdef stats
                Cost before = state.cost();
                auto start = std::chrono::steady_clock::now();
#endif
                if (verbosity) std::cout << "[" << state.cost() << "] round mode=2\n";
                improvement_found = find_omega_improvement(state, k);
#ifdef stats
                ++num_rounds_omega; improvement_omega += state.cost() - before;
                time_omega = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count();
#endif
            }

            // assert(!state.solvable() || bound_graph_is_valid(state, m_edit_state->marked_map(), m_bound_graph));

            rounds_no_improvement = improvement_found ? 0 : rounds_no_improvement + 1;
            // improvement_found || (rounds_no_improvement < 5 && bound_changed)
        } while (state.cost() <= k &&
                 (improvement_found || (rounds_no_improvement < m_max_rounds_no_improvement && bound_changed)));
        if (verbosity) std::cout << "[" << state.cost() << "] end local search\n";
#ifdef stats
        std::cerr << "local_search(state.cost=" << std::setw(4) << cost_before << ", k=" << std::setw(4) << k << ") "
                  << std::setw(4) << cost_before << " -> " << std::setw(4) << state.cost() << "\t "
                  << (state.cost() > k ? "x " : "o ")
                  << "one / two / omega\t "
                  << "# " << num_rounds_one << " " << num_rounds_two << " " << num_rounds_omega << "\t\t"
                  << "+cost " << std::setw(4) << improvement_one << " " << std::setw(4) << improvement_two << " " << std::setw(4) << improvement_omega << "\t\t"
                  << std::setw(8) << time_one << " " << std::setw(8) << time_two << " " << std::setw(8) << time_omega << "\n";
#endif
#undef stats

#ifndef NDEBUG
        if (state.cost() <= k) {
            assert(bound_is_maximal(m_finder, m_edit_state->graph(), m_bound_graph));
        }
#endif
    }

    /**
     * Try to find a (1, 1) swap for the subgraph at the index which improves the lower bound.
     *
     * For the subgraph currently in the lower bound at the given index do the following. List other forbidden subgraphs
     * which are adjacent to the subgraph but not to an other subgraph in the lower bound as candidates.
     * + If a candidate has a larger cost than the original subgraph, choose it instead.
     *
     * Complexity: O((p over 2) * find_near_with_duplicates + |candidates|)
     *
     * @param state
     * @param index
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    bool LocalSearch<SetOfForbiddenSubgraphs>::find_one_improvements(State &state, size_t index) {
        bool found_improvement = false;

        const auto &[subgraph_cost, subgraph] = state.bound(index);
        assert(subgraph_cost == subgraph.calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map()));


        const auto[num_pairs_with_neighbors, num_neighbors_upper_bound] =
            count_neighbors(*m_subgraph_stats, m_edit_state->marked_map(), subgraph);
        if (num_pairs_with_neighbors < 1 || num_neighbors_upper_bound < 1)
            return false; // No improvement possible.

        remove_from_graph(subgraph, m_edit_state->marked_map(), m_bound_graph);


        Cost max_cost = subgraph_cost;
        Subgraph max_subgraph(subgraph);


        const auto pairs = get_pairs(subgraph, m_edit_state->marked_map());

#ifndef NDEBUG
        {
            auto edits = subgraph.non_converting_edits();
            bool touches = std::any_of(edits.begin(), edits.end(), [&](auto uv) {
                return m_bound_graph.has_edge(uv);
            });
            assert(!touches);
        }
#endif

        for (VertexPair uv : pairs) {
            assert(!m_bound_graph.has_edge(uv));

            if (m_subgraph_stats->subgraph_count(uv) > 1) {

                m_finder.find_near(uv, m_edit_state->graph(), m_bound_graph, [&](const Subgraph &neighbor) {
#ifndef NDEBUG
                    {
                        auto edits = neighbor.non_converting_edits();
                        bool touches = std::any_of(edits.begin(), edits.end(), [&](auto xy) {
                            return m_edit_state->graph().has_edge(xy);
                        });
                        assert(!touches);
                    }
#endif
                    Cost n_cost = neighbor.calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map());
                    if (n_cost > max_cost) {
                        found_improvement = true;
                        max_cost = n_cost;
                        max_subgraph = std::move(neighbor);
                    }

                    return subgraph_iterators::IterationControl::Continue;
                });
            }

            // prevent subgraphs including uv to be counted twice
            m_bound_graph.set_edge(uv);
        }

        for (VertexPair uv : pairs)
            m_bound_graph.reset_edge(uv);

        if (max_cost > subgraph_cost && verbosity > 1)
            std::cout << "found (1, 1) swap " << std::setw(4) << max_cost - subgraph_cost << ", " << subgraph_cost
                      << " => " << max_cost << ", " << subgraph << " => " << max_subgraph << "\n";

        insert_into_graph(max_subgraph, m_edit_state->marked_map(), m_bound_graph);
        assert(max_cost == max_subgraph.calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map()));
        state.replace(index, {max_cost, max_subgraph});

        if (max_cost == invalid_cost) {
            state.set_unsolvable();
            std::cerr << "unsolvable one_improvement\n";
        }

        if (state.solvable())
            assert(bound_graph_is_valid(state, m_edit_state->marked_map(), m_bound_graph, m_edit_state->cost_map()));

        return found_improvement;
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
    template<Options::FSG SetOfForbiddenSubgraphs>
    bool LocalSearch<SetOfForbiddenSubgraphs>::find_two_improvement(State &state, size_t index, bool &bound_changed) {
        constexpr Cost invalid_max_cost = std::numeric_limits<Cost>::min();
        constexpr size_t invalid_index = std::numeric_limits<size_t>::max();

        bool found_improvement = false;

        const auto &[subgraph_cost, subgraph] = state.bound(index);
        assert(subgraph_cost == subgraph.calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map()));
        assert(subgraph_cost != invalid_cost);

        // If the subgraph has no neighbors on at least one vertex pair it can be skipped.
        const auto[num_pairs_with_neighbors, num_neighbors_upper_bound] =
            count_neighbors(*m_subgraph_stats, m_edit_state->marked_map(), subgraph);
        if (num_pairs_with_neighbors < 1 || num_neighbors_upper_bound < 1)
            return false; // No improvement possible.

        // Remove the subgraph from m_bound_graph. Either the subgraph or other subgraphs will be reinserted.
        remove_from_graph(subgraph, m_edit_state->marked_map(), m_bound_graph);

        // Candidates are subgraphs which are only adjacent to subgraph but no other subgraph in the lower bound.
        const auto pairs = get_pairs(subgraph, m_edit_state->marked_map());
        auto[candidates, border] = get_candidates(m_finder, pairs, m_edit_state->graph(), m_bound_graph, *m_subgraph_stats);

#ifndef NDEBUG
        {
            auto edits = subgraph.non_converting_edits();
            bool touches = std::any_of(edits.begin(), edits.end(), [&](auto uv) {
                return !m_edit_state->is_marked(uv) && m_bound_graph.has_edge(uv);
            });
            assert(!touches);
        }
#endif

        std::vector<Cost> candidate_costs(candidates.size());
        for (size_t i = 0; i < candidates.size(); ++i)
            candidate_costs[i] = candidates[i].calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map());


        // The information about the best solution.
        Cost max_subgraphs_cost = subgraph_cost;
        std::vector<size_t> max_subgraphs;

        std::vector<size_t> plateau_candidates;

        // for each candidate check if
        //   1. the candidate has a larger cost than the current maximum cost or
        //   2. a pair of two candidates can both be inserted and their cost is larger than the current maximum cost.
        for (size_t pair_i = 0; pair_i < pairs.size(); ++pair_i) {
            for (size_t a_i = border[pair_i]; a_i < border[pair_i + 1]; ++a_i) {
                const auto &a = candidates[a_i];

                insert_into_graph(a, m_edit_state->marked_map(), m_bound_graph);

                Cost max_cost_b = invalid_max_cost;
                size_t max_b_i = invalid_index;

                // for each candidate pair (a, b)
                for (size_t pair_j = pair_i + 1; pair_j < pairs.size(); ++pair_j) {
                    if (m_bound_graph.has_edge(pairs[pair_j])) continue;
                    for (size_t b_i = border[pair_j]; b_i < border[pair_j + 1]; ++b_i) {
                        const auto &b = candidates[b_i];

                        bool inserted = try_insert_into_graph(b, m_edit_state->marked_map(), m_bound_graph);

                        if (inserted) {
                            Cost cost_b = candidate_costs[b_i];
                            if (cost_b > max_cost_b) {
                                max_cost_b = cost_b;
                                max_b_i = b_i;
                            }
                            remove_from_graph(b, m_edit_state->marked_map(), m_bound_graph);
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
                remove_from_graph(a, m_edit_state->marked_map(), m_bound_graph);
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
            insert_into_graph(candidates[a_i], m_edit_state->marked_map(), m_bound_graph);
            assert(candidate_costs[a_i] == candidates[a_i].calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map()));
            state.replace(index, {candidate_costs[a_i], candidates[a_i]});

            if (max_subgraphs.size() == 2) {
                size_t b_i = max_subgraphs[1];

                insert_into_graph(candidates[b_i], m_edit_state->marked_map(), m_bound_graph);
                assert(candidate_costs[b_i] == candidates[b_i].calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map()));
                state.insert({candidate_costs[b_i], candidates[b_i]});
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
                            count_neighbors(*m_subgraph_stats, m_edit_state->marked_map(), candidates[i]);
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

                insert_into_graph(candidates[a_i], m_edit_state->marked_map(), m_bound_graph);
                assert(candidate_costs[a_i] == candidates[a_i].calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map()));
                state.replace(index, {candidate_costs[a_i], candidates[a_i]});
            } else {
                insert_into_graph(subgraph, m_edit_state->marked_map(), m_bound_graph);
            }
        }


        for (size_t pair_i = 0; pair_i < pairs.size(); ++pair_i) {
            if (m_bound_graph.has_edge(pairs[pair_i]))
                continue;
            for (size_t a_i = border[pair_i]; a_i < border[pair_i + 1]; ++a_i) {
                if (candidates[a_i].size() == 0)
                    continue;

                auto inserted = try_insert_into_graph(candidates[a_i], m_edit_state->marked_map(), m_bound_graph);
                if (inserted) {
                    assert(candidate_costs[a_i] == candidates[a_i].calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map()));
                    state.insert({candidate_costs[a_i], candidates[a_i]});
                }
            }
        }


        if (state.solvable()) {
            assert(bound_graph_is_valid(state, m_edit_state->marked_map(), m_bound_graph, m_edit_state->cost_map()));
            assert(bound_is_maximal(m_finder, m_edit_state->graph(), m_bound_graph));
        }

        return found_improvement;
    }

    /**
     * Try to find a (\omega, 1) swap which improves the lower bound.
     *
     * Find all forbidden subgraph in the graph. Check which subgraphs in the lower bound are adjacent to it. If the
     * cost of the single subgraph is larger than the sum of the adjacent subgraphs replace them.
     *
     * Complexity: O(find * (p choose 2))
     *
     * @param state
     * @param k
     * @return
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    bool LocalSearch<SetOfForbiddenSubgraphs>::find_omega_improvement(State &state, Cost k) {
        constexpr size_t invalid_index = std::numeric_limits<size_t>::max();

        assert(state.solvable());

        bool found_improvement = false;

        // VertexPairMap<bool> used(m_bound_graph.size());
        VertexPairMap<size_t> subgraph_index(m_bound_graph.size(), invalid_index);

        // Assign i at all unmarked vertex pairs.
        auto index_assign = [&](const Subgraph &subgraph, size_t i) {
            for (auto uv : subgraph.non_converting_edits()) {
                if (!m_edit_state->is_marked(uv))
                    subgraph_index[uv] = i;
            }
        };

        // Return the indices of the subgraphs already in the lower bound which share a vertex pair with the subgraph.
        auto candidate_neighbors_indices = [&](const Subgraph &subgraph) {
            std::vector<size_t> is;
            for (auto uv : subgraph.non_converting_edits()) {
                is.push_back(subgraph_index[uv]);
            }
            is.erase(std::remove_if(is.begin(), is.end(), [](size_t i) { return i == invalid_index; }), is.end());
            std::sort(is.begin(), is.end());
            is.erase(std::unique(is.begin(), is.end()), is.end());
            return is;
        };

        // populate the subgraph_index map
        for (size_t i = 0; i < state.bound().size(); ++i)
            index_assign(state.bound(i).subgraph, i);

#ifndef NDEBUG
        for (VertexPair xy : Graph::VertexPairs(subgraph_index.size()))
            if (subgraph_index[xy] != invalid_index) {
                auto edits = state.bound(subgraph_index[xy]).subgraph.non_converting_edits();
                auto contains_pair = std::any_of(edits.begin(), edits.end(), [&](auto e) {
                    return e == xy;
                });
                assert(contains_pair);
            }
        for (size_t i = 0; i < state.bound().size(); ++i) {
            for (auto xy : state.bound(i).subgraph.non_converting_edits()) {
                if (!m_edit_state->is_marked(xy))
                    assert(subgraph_index[xy] == i);
            }
        }
#endif

        m_finder.find(m_edit_state->graph(), [&](const Subgraph &subgraph) {

            auto subgraph_cost = subgraph.calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map());

            if (subgraph_cost == invalid_cost) {
                found_improvement = true;
                state.set_unsolvable();
                return subgraph_iterators::IterationControl::Break;
            }
            assert(subgraph_cost != invalid_cost);

            Cost sum = 0;
            size_t count = 0;


            // add the costs of all adjacent subgraphs in the lower bound
            for (auto i : candidate_neighbors_indices(subgraph)) {
                assert(i < state.bound().size());
                sum += state.bound(i).cost;
                count++;
            }

            if (subgraph_cost > sum) {
                // The cost of the single subgraph is larger than the cost of the adjacent subgraphs in the lower bound.
                // Swapping them improves the lower bound.
                found_improvement = true;
                if (verbosity > 1)
                    std::cout << "found (" << count << ", 1) swap " << std::setw(4) << subgraph_cost - sum << ", "
                              << sum << " => " << subgraph_cost << ", ... => " << subgraph << "\n";

                // Remove adjacent subgraphs from the lower bound.
                for (auto uv : subgraph.non_converting_edits()) {
                    auto i = subgraph_index[uv];
                    if (i != invalid_index) {
                        const auto &neighbor = state.bound(i).subgraph;

                        // Overwrite the subgraph_index positions of previously last subgraph.
                        index_assign(state.bound().back().subgraph, i);

                        // Remove neighbor from subgraph_index, m_bound_graph and the bound.
                        index_assign(neighbor, invalid_index);
                        remove_from_graph(neighbor, m_edit_state->marked_map(), m_bound_graph);
                        state.remove(i);
                    }
                }

#ifndef NDEBUG
                if (!bound_graph_is_valid(state, m_edit_state->marked_map(), m_bound_graph, m_edit_state->cost_map())) {
                    for (auto uv : subgraph.non_converting_edits()) {
                        auto i = subgraph_index[uv];
                        if (i != invalid_index) {
                            const auto &neighbor = state.bound(i).subgraph;
                            std::cerr << neighbor << " ";
                        }
                    }
                    std::cerr << "\n";
                    assert(false);
                }
#endif

                // Insert the subgraph into subgraph_index, m_bound_graph and the bound.
                index_assign(subgraph, state.bound().size());
                insert_into_graph(subgraph, m_edit_state->marked_map(), m_bound_graph);
                assert(subgraph_cost == subgraph.calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map()));
                state.insert({subgraph_cost, subgraph});

#ifndef NDEBUG
                for (VertexPair xy : Graph::VertexPairs(subgraph_index.size()))
                    if (subgraph_index[xy] != invalid_index) {
                        auto edits = state.bound(subgraph_index[xy]).subgraph.non_converting_edits();
                        auto contains_pair = std::any_of(edits.begin(), edits.end(), [&](auto e) {
                            return e == xy;
                        });
                        assert(contains_pair);
                    }
                for (size_t j = 0; j < state.bound().size(); ++j) {
                    for (auto xy : state.bound(j).subgraph.non_converting_edits()) {
                        if (!m_edit_state->is_marked(xy))
                            assert(subgraph_index[xy] == j);
                    }
                }
#endif

                // assert(bound_graph_is_valid(state, m_edit_state->marked_map(), m_bound_graph));
            }
            return subgraph_iterators::break_if(state.cost() > k || !state.solvable());
        });

        if (state.cost() <= k && state.solvable()) {
            std::vector<std::pair<Cost, Subgraph>> remaining_subgraphs;
            m_finder.find(m_edit_state->graph(), m_bound_graph, [&](const Subgraph &subgraph) {
                auto cost = subgraph.calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map());
                remaining_subgraphs.emplace_back(cost, subgraph);
                return subgraph_iterators::IterationControl::Continue;
            });

            for (auto &&[cost, subgraph] : remaining_subgraphs) {
                bool inserted = try_insert_into_graph(subgraph, m_edit_state->marked_map(), m_bound_graph);
                if (inserted) {
                    assert(cost == subgraph.calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map()));
                    state.insert({cost, subgraph});
                }
            }
        }

        if (state.cost() <= k && state.solvable()) {
            assert(bound_graph_is_valid(state, m_edit_state->marked_map(), m_bound_graph, m_edit_state->cost_map()));
            assert(bound_is_maximal(m_finder, m_edit_state->graph(), m_bound_graph));
        }

        return found_improvement;
    }

    /**
     * Insert unmarked vertex pairs of "subgraph" into "graph".
     *
     * @param subgraph
     * @param forbidden
     * @param graph
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    void LocalSearch<SetOfForbiddenSubgraphs>::insert_into_graph(const LocalSearch::Subgraph &subgraph,
                                                                 const VertexPairMap<bool> &marked, Graph &graph) {
        for (auto uv : subgraph.non_converting_edits()) {
            if (!marked[uv]) {
                assert(!graph.has_edge(uv));
                graph.set_edge(uv);
            }
        }
    }


    /**
     * Remove unmarked vertex pairs of "subgraph" from "graph".
     *
     * @param subgraph
     * @param forbidden
     * @param graph
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    void LocalSearch<SetOfForbiddenSubgraphs>::remove_from_graph(const Subgraph &subgraph,
                                                                 const VertexPairMap<bool> &marked, Graph &graph) {
        for (auto uv : subgraph.non_converting_edits()) {
            if (!marked[uv]) {
                assert(graph.has_edge(uv));
                graph.reset_edge(uv);
            }
        }
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
    template<Options::FSG SetOfForbiddenSubgraphs>
    bool LocalSearch<SetOfForbiddenSubgraphs>::try_insert_into_graph(const Subgraph &subgraph,
                                                                     const VertexPairMap<bool> &marked, Graph &graph) {
        auto edits = subgraph.non_converting_edits();
        bool touches = std::any_of(edits.begin(), edits.end(), [&](auto uv) {
            return !marked[uv] && graph.has_edge(uv);
        });


        if (!touches) {
            for (auto uv : subgraph.non_converting_edits()) {
                if (!marked[uv])
                    graph.set_edge(uv);
            }
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
    template<Options::FSG SetOfForbiddenSubgraphs>
    std::pair<std::vector<typename LocalSearch<SetOfForbiddenSubgraphs>::Subgraph>, std::vector<size_t>>
    LocalSearch<SetOfForbiddenSubgraphs>::get_candidates(Finder &finder, const std::vector<VertexPair> &pairs,
                                                         const Graph &graph,
                                                         Graph &bound_graph,
                                                         const SubgraphStats<SetOfForbiddenSubgraphs> &subgraph_stats) {

        // Precondition: subgraph is removed from bound_graph.
#ifndef NDEBUG
        for (auto uv : pairs) {
            assert(!bound_graph.has_edge(uv));
        }
#endif

        std::vector<Subgraph> candidates;
        std::vector<size_t> border(pairs.size() + 1);

        for (size_t i = 0; i < pairs.size(); ++i) {
            VertexPair uv = pairs[i];
            assert(!bound_graph.has_edge(uv));

            if (subgraph_stats.subgraph_count(uv) > 1) {
                finder.find_near(uv, graph, bound_graph, [&](const Subgraph &neighbor) {
#ifndef NDEBUG
                    for (auto xy : neighbor.non_converting_edits()) {
                        assert(!bound_graph.has_edge(xy));
                    }
#endif
                    candidates.push_back(neighbor);
                    return subgraph_iterators::IterationControl::Continue;
                });
            }
            border[i + 1] = candidates.size();

            // prevent subgraphs including uv to be counted twice
            bound_graph.set_edge(uv);
        }

        // reset bound_graph
        for (VertexPair uv : pairs) {
            assert(bound_graph.has_edge(uv));
            bound_graph.reset_edge(uv);
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
    template<Options::FSG SetOfForbiddenSubgraphs>
    std::vector<VertexPair> LocalSearch<SetOfForbiddenSubgraphs>::get_pairs(const Subgraph &subgraph,
                                                                            const VertexPairMap<bool> &marked) {
        std::vector<VertexPair> pairs;
        for (auto uv : subgraph.non_converting_edits()) {
            if (!marked[uv])
                pairs.push_back(uv);
        }
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
    template<Options::FSG SetOfForbiddenSubgraphs>
    std::tuple<size_t, size_t> LocalSearch<SetOfForbiddenSubgraphs>::count_neighbors(
            const SubgraphStats<SetOfForbiddenSubgraphs> &subgraph_stats, const VertexPairMap<bool> &marked,
            const Subgraph &subgraph) {
        size_t num_pairs = 0, num_neighbors_ub = 0;
        for (auto uv : subgraph.non_converting_edits()) {
            if (!marked[uv]) {
                size_t nn = subgraph_stats.subgraph_count(uv) - 1;
                num_neighbors_ub += nn;
                if (nn > 0) ++num_pairs;
            }
        }
        return {num_pairs, num_neighbors_ub};
    }

    template
    class LocalSearch<Options::FSG::C4P4>;

}
