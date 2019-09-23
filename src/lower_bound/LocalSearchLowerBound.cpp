//
// Created by jonas on 25.07.19.
//


#include "LocalSearchLowerBound.h"

#include <chrono>


/**
 * Returns a lower bound on the editing cost. The state keeps a maximal lower bound. Local search tries to improve
 * bound by performing (1, 1), (1, 2) or (\omega, 1) swaps on the forbidden subgraphs to improve the lower bound.
 *
 * @param k Remaining editing cost.
 * @return A lower bound on the costs required to solve the current instance.
 */
Cost LocalSearchLowerBound::calculate_lower_bound(Cost k) {
    auto &state = current_state();
    if (!state.solvable()) return state.cost();

    // The bound is beeing kept maximal. When the cost is zero, no improvement can be made.
    if (state.cost() == 0) return 0;

    state.recalculate(m_marked, m_costs);
    state.initialize_bound_graph(m_marked, m_bound_graph);
    state.assert_valid(m_marked, m_bound_graph, m_costs);

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
void LocalSearchLowerBound::initialize(Cost /*k*/) {
    using Element = State::Element;
    m_states.push_back(std::make_unique<State>());
    State &state = *m_states.back();

    m_bound_graph.clearEdges();


    std::vector<Element> subgraphs;
    finder->find([&](Subgraph &&subgraph) {

        Cost min_cost = get_subgraph_cost(subgraph, m_marked, m_costs);
        subgraphs.push_back({min_cost, std::move(subgraph)});

        // If a subgraph exists with fully marked vertex pairs, the instance is not solvable.
        return min_cost == invalid_cost;
    });

    std::sort(subgraphs.begin(), subgraphs.end(),
              [](const Element &lhs, const Element &rhs) { return lhs.cost > rhs.cost; });

    if (!subgraphs.empty() && subgraphs[0].cost == invalid_cost) {
        state.set_unsolvable();
        return;
    }

    for (auto&[cost, subgraph] : subgraphs) {
        bool inserted = try_insert_into_graph(subgraph, m_marked, m_bound_graph);

        if (inserted) {
            assert(cost != invalid_cost);
            state.insert({cost, std::move(subgraph)});
        }
    }
}

void LocalSearchLowerBound::remove_subgraphs_from_bound(State& state, VertexPair uv) {
    assert(state.solvable());
    // state.assert_valid(m_marked, m_bound_graph, m_costs);

    for (size_t i = 0; i < state.bound().size();) {
        const auto &[cost, subgraph] = state.bound(i);
        const auto vertices = subgraph.vertices();

        bool has_u = std::any_of(vertices.begin(), vertices.end(), [&](Vertex x) { return x == uv.u; });
        bool has_v = std::any_of(vertices.begin(), vertices.end(), [&](Vertex x) { return x == uv.v; });

        if (has_u && has_v) {
            // std::cout << "before_mark_and_edit remove " << subgraph;
            state.remove(i);
            // std::cout << "\t => " << state << "\n";
        } else {
            ++i;
        }
    }
}

void LocalSearchLowerBound::insert_subgraphs_into_bound(State& state, VertexPair uv) {
    assert(state.solvable());

    // TODO: Check if needed
    state.initialize_bound_graph(m_marked, m_bound_graph);

    std::vector<std::pair<Cost, Subgraph>> subgraphs;

    // The finder iterates over subgraphs having u and v as vertices.
    finder->find_near(uv, m_bound_graph, [&](Subgraph &&subgraph) {
        Cost min_cost = get_subgraph_cost(subgraph, m_marked, m_costs);
        subgraphs.emplace_back(min_cost, std::move(subgraph));
        return false;
    });

    std::sort(subgraphs.begin(), subgraphs.end(),
              [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

    // if one forbidden subgraph has no unmarked pairs, the instance is unsolvable.
    if (!subgraphs.empty() && subgraphs[0].first == invalid_cost) {
        state.set_unsolvable();
        return;
    }

    for (auto &[cost, subgraph] : subgraphs) {

        // In the for loop a subgraph can be invalidated when bound_graph is modified.
        // Therefore it has to checked whether subgraph is still valid and does not share a vertex pair with bound_graph.
        const auto vertexPairs = subgraph.vertexPairs();
        bool touches_bound = std::any_of(vertexPairs.begin(), vertexPairs.end(), [&](VertexPair xy) {
            return !m_marked[xy] && m_bound_graph.hasEdge(xy);
        });

        if (!touches_bound) {
            // std::cout << "after_mark_and_edit insert " << subgraph;
            insert_into_graph(subgraph, m_marked, m_bound_graph);
            state.insert({cost, std::move(subgraph)});
            // std::cout << "\t => " << state << "\n";
        }
    }
}

/**
 * Removes subgraphs from state.bound which have the vertex pair uv.
 *
 * @param uv
 */
void LocalSearchLowerBound::before_mark_and_edit(VertexPair uv) {
    auto &state = current_state();
    remove_subgraphs_from_bound(state, uv);
}

void LocalSearchLowerBound::before_mark(VertexPair uv) {
    auto &state = parent_state();
    remove_subgraphs_from_bound(state, uv);
}

/**
 * Find forbidden subgraphs near uv and insert them. Subgraphs with higher minimum editing cost are preferred.
 *
 * @param uv
 */
void LocalSearchLowerBound::after_mark_and_edit(VertexPair uv) {
    auto &state = current_state();
    insert_subgraphs_into_bound(state, uv);
}

/**
 * Tries to optimise the bound. Stops if the lower bound is larger than k or the bound does not improve in five consecutive rounds.
 *
 * @param state
 * @param k
 */
void LocalSearchLowerBound::local_search(State &state, Cost k) {
    constexpr static bool use_one_improvement = false;
    constexpr static bool use_two_improvement = true;
    constexpr static bool use_omega_improvement = true;
#define stats
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

        rounds_no_improvement = improvement_found ? 0 : rounds_no_improvement + 1;
        // improvement_found || (rounds_no_improvement < 5 && bound_changed)
    } while (state.cost() <= k && (improvement_found || (rounds_no_improvement < m_max_rounds_no_improvement && bound_changed)));
    if (verbosity) std::cout << "[" << state.cost() << "] end local search\n";
#ifdef stats
    std::cerr << "local_search(state.cost=" << std::setw(4) << cost_before << ", k=" << std::setw(4) << k << ") "
              << std::setw(4) << cost_before << " -> " << std::setw(4) << state.cost() << "\t "
              << "one / two / omega\t "
              << "# " << num_rounds_one << " " << num_rounds_two << " " << num_rounds_omega << "\t\t"
              << "+cost " << std::setw(4) << improvement_one << " " << std::setw(4) << improvement_two << " " << std::setw(4) << improvement_omega << "\t\t"
              << std::setw(6) << time_one << "\t " << std::setw(6) << time_two << "\t " << std::setw(6) << time_omega << "\n";
#endif
#undef stats
}

/**
 * Try to find a (1, 1) swap for the subgraph at the index which improves the lower bound.
 *
 * For the subgraph currently in the lower bound at the given index do the following. List other forbidden subgraphs
 * which are adjacent to the subgraph but not to an other subgraph in the lower bound as candidates.
 * + If a candidate has a larger cost than the original subgraph, choose it instead.
 *
 * Complexity: O((p over 2) * find_near + |candidates|)
 *
 * @param state
 * @param index
 */
bool LocalSearchLowerBound::find_one_improvements(State &state, size_t index) {
    bool found_improvement = false;

    const auto &[subgraph_cost, subgraph] = state.bound(index);
    assert(subgraph_cost == get_subgraph_cost(subgraph, m_marked, m_costs));


    const auto [num_pairs_with_neighbors, num_neighbors_upper_bound] = count_neighbors(m_subgraph_stats, m_marked, subgraph);
    if (num_pairs_with_neighbors < 1 || num_neighbors_upper_bound < 1)
        return false; // No improvement possible.

    remove_from_graph(subgraph, m_marked, m_bound_graph);


    Cost max_cost = subgraph_cost;
    Subgraph max_subgraph(subgraph);


    const auto pairs = get_pairs(subgraph, m_marked);

#ifndef NDEBUG
    {
        bool touches = false;
        for (VertexPair uv : subgraph.vertexPairs())
            if (m_bound_graph.hasEdge(uv))
                touches = true;
        assert(!touches);
    }
#endif

    for (VertexPair uv : pairs) {
        assert(!m_bound_graph.hasEdge(uv));

        if (m_subgraph_stats.subgraphCount(uv) > 1) {

            finder->find_near(uv, m_bound_graph, [&](Subgraph &&neighbor) {
#ifndef NDEBUG
                {
                    auto vp = neighbor.vertexPairs();
                    assert(std::none_of(vp.begin(), vp.end(),
                                        [&](VertexPair xy) { return m_bound_graph.hasEdge(xy); }));
                }
#endif
                Cost n_cost = get_subgraph_cost(neighbor, m_marked, m_costs);
                if (n_cost > max_cost) {
                    found_improvement = true;
                    max_cost = n_cost;
                    max_subgraph = std::move(neighbor);
                }

                return false;
            });
        }

        // prevent subgraphs including uv to be counted twice
        m_bound_graph.setEdge(uv);
    }

    for (VertexPair uv : pairs)
        m_bound_graph.clearEdge(uv);

    if (max_cost > subgraph_cost && verbosity > 1)
        std::cout << "found (1, 1) swap " << std::setw(4) << max_cost - subgraph_cost << ", " << subgraph_cost
                  << " => " << max_cost << ", " << subgraph << " => " << max_subgraph << "\n";

    insert_into_graph(max_subgraph, m_marked, m_bound_graph);
    state.replace(index, {max_cost, std::move(max_subgraph)});

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
 * Complexity: O((p over 2) * find_near + |candidates|^2)
 *
 * @param state
 * @param index
 * @param gen
 * @param bound_changed
 * @return
 */
bool LocalSearchLowerBound::find_two_improvement(State &state, size_t index, bool &bound_changed) {
    constexpr Cost invalid_max_cost = std::numeric_limits<Cost>::min();
    constexpr size_t invalid_index = std::numeric_limits<size_t>::max();

    bool found_improvement = false;

    const auto &[subgraph_cost, subgraph] = state.bound(index);
    assert(subgraph_cost == get_subgraph_cost(subgraph, m_marked, m_costs));
    assert(subgraph_cost != invalid_cost);

    // If the subgraph has no neighbors on at least one vertex pair it can be skipped.
    const auto [num_pairs_with_neighbors, num_neighbors_upper_bound] = count_neighbors(m_subgraph_stats, m_marked, subgraph);
    if (num_pairs_with_neighbors < 1 || num_neighbors_upper_bound < 1)
        return false; // No improvement possible.

    // Remove the subgraph from m_bound_graph. Either the subgraph or other subgraphs will be reinserted.
    remove_from_graph(subgraph, m_marked, m_bound_graph);

    // Candidates are subgraphs which are only adjacent to subgraph but no other subgraph in the lower bound.
    const auto pairs = get_pairs(subgraph, m_marked);
    auto[candidates, border] = get_candidates(*finder, pairs, m_bound_graph, m_subgraph_stats);

#ifndef NDEBUG
    {
        bool touches = false;
        for (VertexPair uv : subgraph.vertexPairs())
            if (!m_marked[uv] && m_bound_graph.hasEdge(uv)) touches = true;
        assert(!touches);
    }
#endif

    std::vector<Cost> candidate_costs(candidates.size());
    for (size_t i = 0; i < candidates.size(); ++i)
        candidate_costs[i] = get_subgraph_cost(candidates[i], m_marked, m_costs);


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

            insert_into_graph(a, m_marked, m_bound_graph);

            Cost max_cost_b = invalid_max_cost;
            size_t max_b_i = invalid_index;

            // for each candidate pair (a, b)
            for (size_t pair_j = pair_i + 1; pair_j < pairs.size(); ++pair_j) {
                if (m_bound_graph.hasEdge(pairs[pair_j])) continue;
                for (size_t b_i = border[pair_j]; b_i < border[pair_j + 1]; ++b_i) {
                    const Subgraph &b = candidates[b_i];

                    bool inserted = try_insert_into_graph(b, m_marked, m_bound_graph);

                    if (inserted) {
                        Cost cost_b = candidate_costs[b_i];
                        if (cost_b > max_cost_b) {
                            max_cost_b = cost_b;
                            max_b_i = b_i;
                        }
                        remove_from_graph(b, m_marked, m_bound_graph);
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
            remove_from_graph(a, m_marked, m_bound_graph);
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
        insert_into_graph(candidates[a_i], m_marked, m_bound_graph);
        state.replace(index, {candidate_costs[a_i], std::move(candidates[a_i])});

        if (max_subgraphs.size() == 2) {
            size_t b_i = max_subgraphs[1];

            insert_into_graph(candidates[b_i], m_marked, m_bound_graph);
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
                    auto [num_pairs, num_covered_ub] = count_neighbors(m_subgraph_stats, m_marked, candidates[i]);
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

            insert_into_graph(candidates[a_i], m_marked, m_bound_graph);
            state.replace(index, {candidate_costs[a_i], std::move(candidates[a_i])});
        } else {
            insert_into_graph(subgraph, m_marked, m_bound_graph);
        }
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
bool LocalSearchLowerBound::find_omega_improvement(State &state, Cost k) {
    constexpr size_t invalid_index = std::numeric_limits<size_t>::max();

    bool found_improvement = false;

    // VertexPairMap<bool> used(m_bound_graph.size());
    VertexPairMap<size_t> subgraph_index(m_bound_graph.size(), invalid_index);

    // Assign i at all unmarked vertex pairs.
    auto index_assign = [&](const Subgraph &subgraph, size_t i) {
        for (VertexPair uv : subgraph.vertexPairs())
            if (!m_marked[uv])
                subgraph_index[uv] = i;
    };

    // Return the indices of the subgraphs already in the lower bound which share a vertex pair with the subgraph.
    auto candidate_neighbors_indices = [&](const Subgraph &subgraph) {
        std::vector<size_t> is;
        for (VertexPair uv : subgraph.vertexPairs())
            is.push_back(subgraph_index[uv]);
        is.erase(std::remove_if(is.begin(), is.end(), [](size_t i) { return i == invalid_index; }), is.end());
        std::sort(is.begin(), is.end());
        is.erase(std::unique(is.begin(), is.end()), is.end());
        return is;
    };

    // populate the subgraph_index map
    for (size_t i = 0; i < state.bound().size(); ++i)
        index_assign(state.bound(i).subgraph, i);


    finder->find([&](Subgraph &&subgraph) {
        auto subgraph_cost = get_subgraph_cost(subgraph, m_marked, m_costs);

        if (subgraph_cost == invalid_cost) {
            state.set_unsolvable();
            return true;
        }
        assert(subgraph_cost != invalid_cost);

        Cost sum = 0;
        size_t count = 0;


        // add the costs of all adjacent subgraphs in the lower bound
        for (auto i : candidate_neighbors_indices(subgraph)) {
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
            for (VertexPair uv : subgraph.vertexPairs()) {
                auto i = subgraph_index[uv];
                if (i != invalid_index) {
                    const auto &neighbor = state.bound(i).subgraph;

                    // Remove neighbor from subgraph_index, m_bound_graph and the bound.
                    index_assign(neighbor, invalid_index);
                    remove_from_graph(neighbor, m_marked, m_bound_graph);
                    state.remove(i);

                    // Overwrite the subgraph_index positions of previously last subgraph.
                    index_assign(state.bound(i).subgraph, i);
                }
            }

            // Insert the subgraph into subgraph_index, m_bound_graph and the bound.
            index_assign(subgraph, state.bound().size());
            insert_into_graph(subgraph, m_marked, m_bound_graph);
            state.insert({subgraph_cost, std::move(subgraph)});
        }
        return state.cost() > k;
    });

    return found_improvement;
}

/**
 * Insert unmarked vertex pairs of "subgraph" into "graph".
 *
 * @param subgraph
 * @param forbidden
 * @param graph
 */
void LocalSearchLowerBound::insert_into_graph(const Subgraph &subgraph, const VertexPairMap<bool> &marked, Graph &graph) {
    for (VertexPair uv : subgraph.vertexPairs()) {
        if (!marked[uv]) {
            assert(!graph.hasEdge(uv));
            graph.setEdge(uv);
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
void LocalSearchLowerBound::remove_from_graph(const Subgraph &subgraph, const VertexPairMap<bool> &marked, Graph &graph) {
    for (VertexPair uv : subgraph.vertexPairs()) {
        if (!marked[uv]) {
            assert(graph.hasEdge(uv));
            graph.clearEdge(uv);
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
bool LocalSearchLowerBound::try_insert_into_graph(const Subgraph &subgraph, const VertexPairMap<bool> &marked, Graph &graph) {
    bool touches = false;
    for (VertexPair uv : subgraph.vertexPairs())
        if (!marked[uv] && graph.hasEdge(uv)) {
            touches = true;
            break;
        }


    if (!touches) {
        for (VertexPair uv : subgraph.vertexPairs())
            if (!marked[uv])
                graph.setEdge(uv);

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
std::pair<std::vector<Subgraph>, std::vector<size_t>> LocalSearchLowerBound::get_candidates(FinderI &finder,
        const std::vector<VertexPair> &pairs, Graph &bound_graph, const SubgraphStats &subgraph_stats) {

    // Precondition: subgraph is removed from bound_graph.
    assert(std::none_of(pairs.begin(), pairs.end(), [&](VertexPair uv) { return bound_graph.hasEdge(uv); }));

    std::vector<Subgraph> candidates;
    std::vector<size_t> border(pairs.size() + 1);

    for (size_t i = 0; i < pairs.size(); ++i) {
        VertexPair uv = pairs[i];
        assert(!bound_graph.hasEdge(uv));

        if (subgraph_stats.subgraphCount(uv) > 1) {
            finder.find_near(uv, bound_graph, [&](Subgraph &&neighbor) {
#ifndef NDEBUG
                auto vp = neighbor.vertexPairs();
                assert(std::none_of(vp.begin(), vp.end(), [&](VertexPair xy) { return bound_graph.hasEdge(xy); }));
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
std::vector<VertexPair> LocalSearchLowerBound::get_pairs(const Subgraph &subgraph, const VertexPairMap<bool> &marked) {
    std::vector<VertexPair> pairs;
    for (VertexPair uv : subgraph.vertexPairs())
        if (!marked[uv])
            pairs.push_back(uv);
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
std::tuple<size_t, size_t> LocalSearchLowerBound::count_neighbors(const SubgraphStats& subgraph_stats,
                                                  const VertexPairMap<bool> &marked, const Subgraph& subgraph) {
    size_t num_pairs = 0, num_neighbors_ub = 0;
    for (VertexPair uv : subgraph.vertexPairs()) {
        if (!marked[uv]) {
            size_t nn = subgraph_stats.subgraphCount(uv) - 1;
            num_neighbors_ub += nn;
            if (nn > 0) ++num_pairs;
        }
    }
    return {num_pairs, num_neighbors_ub};
}
