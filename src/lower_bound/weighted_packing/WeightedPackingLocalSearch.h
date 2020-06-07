#ifndef WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKINGLOCALSEARCH_H
#define WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKINGLOCALSEARCH_H

#include "WeightedPacking.h"

#include "robin_hood.h"
#include <unordered_set>


class WeightedPackingLocalSearch : public LowerBoundI {
    class State {
    private:
        robin_hood::unordered_map<Subgraph, Cost> m_subgraphs_in_packing;
        bool m_solvable = true;
    public:
        [[nodiscard]] auto &subgraphs() {
            return m_subgraphs_in_packing;
        }

        [[nodiscard]] const auto &subgraphs() const {
            return m_subgraphs_in_packing;
        }

        void set_unsolvable() {
            m_solvable = false;
        }

        void set_solvable() {
            m_solvable = true;
        }

        [[nodiscard]] bool is_solvable() const {
            return m_solvable;
        }

#ifndef NDEBUG

        void packing_is_same_as_initialized_by_state(const WeightedPacking &packing) {
            WeightedPacking other(packing);
            other.clear();

            for (const auto &[subgraph, cost] : m_subgraphs_in_packing) {
                other.subtract_from_potential(subgraph, cost);
            }

            assert(other.cost() == packing.cost());
            for (auto uv : other.depleted_graph().vertexPairs()) {
                if (other.potential(uv) != packing.potential(uv)) std::cerr << uv << std::endl;
                assert(other.potential(uv) == packing.potential(uv));
                if (other.is_depleted(uv) != packing.is_depleted(uv)) std::cerr << uv << std::endl;
                assert(other.is_depleted(uv) == packing.is_depleted(uv));
            }
        }

#endif

        void clear() {
            m_subgraphs_in_packing.clear();
            m_solvable = true;
        }

        friend std::ostream &operator<<(std::ostream &os, const State &state) {
            for (const auto&[k, v] : state.subgraphs()) {
                std::cout << k << ": " << v << ", ";
            }
            return os << "\n";
        }
    };

    const Graph &m_graph;
    const VertexPairMap<Cost> &m_costs;
    const VertexPairMap<bool> &m_marked;

    WeightedPacking m_packing;
    std::vector<State> m_states;

    std::mt19937_64 m_gen;

    int verbosity = 0;

public:
    WeightedPackingLocalSearch(const Instance &instance, const VertexPairMap<bool> &marked,
                               const SubgraphStats &subgraph_stats, std::shared_ptr<FinderI> finder_ref) :
            LowerBoundI(finder_ref), m_graph(instance.graph), m_costs(instance.costs), m_marked(marked),
            m_packing(instance, marked, subgraph_stats, std::move(finder_ref)), m_gen(0) {
        m_states.emplace_back();
    }

    void push_state(Cost /*k*/) override {
        if (verbosity > 0) std::cout << "push_state(...) " << m_states.size() << "\n";
        assert(!m_states.empty());
        m_states.emplace_back(m_states.back());
    }

    void pop_state() override {
        if (verbosity > 0) std::cout << "pop_state() " << m_states.size() << "\n";
        assert(!m_states.empty());
        m_states.pop_back();
        // initialize_bound_graph(*finder, current_state(), m_marked, m_bound_graph);
        auto &state = current_state();
        m_packing.clear();
        for (const auto &[subgraph, cost] : state.subgraphs()) {
            m_packing.subtract_from_potential(subgraph, cost);
        }
    }

    State &current_state() {
        assert(!m_states.empty());
        return m_states.back();
    }

    State &parent_state() {
        assert(m_states.size() > 1);
        return m_states[m_states.size() - 2];
    }

    /**
     * Ensure that the state remains maximal by increasing costs of incident subgraphs in the packing.
     *
     * @param uv
     */
    void after_mark(VertexPair uv) override {
        assert(m_marked[uv]);
        if (verbosity > 0) std::cout << "after_mark(" << uv << ")\n";
        auto &state = current_state();

        assert(m_packing.is_valid());
        if (state.is_solvable()) {
            assert(m_packing.is_maximal());
        }

        m_packing.restore_potential(uv);

        if (!state.is_solvable()) return;

#ifndef NDEBUG
        state.packing_is_same_as_initialized_by_state(m_packing);
#endif

        // TODO: Early exit if uv was not depleted.

        // III
        // Note: Can be optimized by inserting subgraphs in a greedy order.
        bool unsolvable = finder->find_near_with_duplicates(uv, m_graph, m_packing.depleted_graph(),
                [&](Subgraph &&subgraph) {
            assert(subgraph.contains(uv));
            auto cost = m_packing.calculate_min_cost(subgraph);
            if (cost == 0) return false;
            if (cost == invalid_cost) {
                return true;
            }
            m_packing.subtract_from_potential(subgraph, cost);
            state.subgraphs()[std::move(subgraph)] += cost;
            return false;
        });

        if (unsolvable) {
            state.set_unsolvable();
            return;
        }

        assert(m_packing.is_maximal());
    }


    /**
     * Replace destroyed subgraphs by newly created ones.
     * @param uv
     */
    void after_edit(VertexPair uv) override {
        assert(m_marked[uv]);
        if (verbosity > 0) std::cout << "after_edit(" << uv << ")\n";
        auto &state = current_state();

        if (verbosity > 0) std::cout << "solvable = " << state.is_solvable() << "\n";
        //if (!state.is_solvable() && state.was_made_solvable_by(uv))
        state.set_solvable();

        assert(m_packing.is_valid());
#ifndef NDEBUG
        state.packing_is_same_as_initialized_by_state(m_packing);
#endif

        assert(m_packing.potential(uv) == m_costs[uv]);

        // Remove subgraphs which are destroyed.
        // Note: This can "free" subgraphs which now can be inserted. These are at pairs which were previously depleted.
        std::vector<VertexPair> pairs;
        if (!m_packing.is_depleted(uv))
            pairs.push_back(uv);

        for (auto it = state.subgraphs().begin(); it != state.subgraphs().end();) {
            const auto &[subgraph, cost] = *it;
            if (finder->for_all_conversionless_edits(subgraph, [&](auto xy) { return xy == uv; })) {
                finder->for_all_conversionless_edits(subgraph, [&](auto xy) {
                    if (!m_marked[xy] && m_packing.is_depleted(xy)) {
                        pairs.push_back(xy);
                    }
                    return false;
                });
                if (verbosity > 0) std::cout << "removing " << subgraph << " " << cost << "\n";
                m_packing.add_to_potential(subgraph, cost);
                it = state.subgraphs().erase(it);
            } else {
                ++it;
            }
        }
        if (verbosity > 0) {
            if (!pairs.empty())
                std::cout << "### pairs :";
            for (auto xy : pairs)
                std::cout << " " << xy;
            std::cout << std::endl;
        }

#ifndef NDEBUG
        for (auto xy : pairs) {
            assert(!m_packing.is_depleted(xy));
        }
#endif

        // Note: Can be optimized by inserting subgraphs in a greedy order.
        auto subgraphs = m_packing.get_incident_subgraphs(pairs);
        for (auto &&subgraph : subgraphs) {
            auto cost = m_packing.calculate_min_cost(subgraph);
            if (cost == 0) continue;
            if (cost == invalid_cost) {
                state.set_unsolvable();
                return;
            }
            if (verbosity > 0) std::cout << "inserting " << subgraph << " " << cost << std::endl;
            m_packing.subtract_from_potential(subgraph, cost);
            state.subgraphs()[std::move(subgraph)] += cost;
        }

        assert(m_packing.is_maximal());
        assert(m_packing.is_valid());
#ifndef NDEBUG
        state.packing_is_same_as_initialized_by_state(m_packing);
#endif
    }

    void initialize(Cost k) override {
        if (verbosity > 0) std::cout << "initialize\n";
        auto &state = current_state();
        greedy_initialize(state, k);
#ifndef NDEBUG
        state.packing_is_same_as_initialized_by_state(m_packing);
#endif
    }

    Cost calculate_lower_bound(Cost k) override {
        if (verbosity > 0) std::cout << "calculate_lower_bound\n";
        auto &state = current_state();
        if (state.is_solvable()) {
            assert(m_packing.is_maximal());
            assert(m_packing.is_valid());
#ifndef NDEBUG
            state.packing_is_same_as_initialized_by_state(m_packing);
#endif
            local_search(state, k);
            return m_packing.cost();
        } else {
            return invalid_cost;
        }
    }

    Cost calculate_lower_bound_no_edit_branch() override {
        if (verbosity > 0) std::cout << "calculate_lower_bound_no_edit_branch\n";
        auto &state = current_state();
        if (state.is_solvable()) {
            assert(m_packing.is_maximal());
            assert(m_packing.is_valid());
#ifndef NDEBUG
            state.packing_is_same_as_initialized_by_state(m_packing);
#endif
            return m_packing.cost();
        } else {
            return invalid_cost;
        }
    }

    void local_search(State &state, Cost k) {
        if (verbosity > 0) std::cout << "local_search\n";
        size_t max_iter = 0;
        size_t max_rounds_no_improvement = 5;
        size_t num_rounds_no_improvement = 0;
        std::vector<std::pair<Subgraph, Cost>> removed_subgraphs, inserted_subgraphs;

        for (size_t iter = 0; iter < max_iter; ++iter) {
            bool changed_in_round = false, improved_in_round = false;
            removed_subgraphs.clear();
            inserted_subgraphs.clear();

            for (const auto &[subgraph, cost] : state.subgraphs()) {
                auto[changed, improved] = \
                    find_one_two_improvement(subgraph, cost, removed_subgraphs, inserted_subgraphs);
                changed_in_round |= changed;
                improved_in_round |= improved;
            }

            // Update state. During a round, the state is stale and updates are only represented in m_packing.
            // Insert before remove, because both could cancel each other.
            for (auto&&[subgraph, cost] : inserted_subgraphs) {
                // If the subgraph is not yet in the map, the initial value is 0.
                state.subgraphs()[std::move(subgraph)] += cost;
            }
            for (auto&&[subgraph, cost] : removed_subgraphs) {
                // If the subgraph is fully removed, i.e. its cost is zero, it is removed from the map.
                auto it = state.subgraphs().find(subgraph);
                assert(it->second >= cost);
                it->second -= cost;
                if (it->second == 0) {
                    state.subgraphs().erase(it);
                }
            }

            num_rounds_no_improvement = improved_in_round ? 0 : num_rounds_no_improvement + 1;
            if (!changed_in_round ||
                num_rounds_no_improvement > max_rounds_no_improvement ||
                m_packing.cost() > k) {
                break;
            }
        }
    }

    void greedy_initialize(State &state, Cost k) {
        if (verbosity > 0) std::cout << "greedy_initialize\n";
        std::vector<std::pair<Subgraph, Cost>> subgraph_heap; // Note: Could turn this into an attribute and call clear() at beginning. This would prevent reallocating memory.

        auto comp = [](const auto &lhs, const auto &rhs) { return lhs.second < rhs.second; };

        bool unsolvable = finder->find_with_duplicates(m_graph, [&](Subgraph &&subgraph) {
            Cost initial_min_cost = m_packing.calculate_min_cost(subgraph);

            subgraph_heap.emplace_back(std::move(subgraph), initial_min_cost);

            return initial_min_cost > k;
        });

        if (unsolvable) {
            state.set_unsolvable();
            return;
        }

        std::make_heap(subgraph_heap.begin(), subgraph_heap.end(), comp);

        while (!subgraph_heap.empty() && m_packing.cost() <= k) {
            std::pop_heap(subgraph_heap.begin(), subgraph_heap.end(), comp);
            auto &[subgraph, cost] = subgraph_heap.back();

            Cost current_min_cost = m_packing.calculate_min_cost(subgraph);

            if (current_min_cost == 0) {
                // The subgraph cannot be inserted.
                subgraph_heap.pop_back();
            } else if (current_min_cost == cost) {
                // The editing cost is maximal from all remaining subgraphs in the queue.
                // Update packing and state.
                m_packing.subtract_from_potential(subgraph, cost);
                // state.subgraphs().emplace(std::move(subgraph_heap.back())); // Alternative to next line.
                state.subgraphs()[std::move(subgraph)] = cost;

                subgraph_heap.pop_back();
            } else {
                // The editing cost is no longer up to date. The subgraph may be inserted in the future.
                subgraph_heap.back().second = current_min_cost;
                std::push_heap(subgraph_heap.begin(), subgraph_heap.end(), comp);
            }
        }

        if (m_packing.cost() <= k) {
            assert(m_packing.is_maximal());
        }
    }


    std::tuple<bool, bool> find_one_two_improvement(const Subgraph &x, Cost x_cost,
                                                    std::vector<std::pair<Subgraph, Cost>> &removed_subgraphs,
                                                    std::vector<std::pair<Subgraph, Cost>> &inserted_subgraphs) {
        if (verbosity > 0) std::cout << "find_one_two_improvement " << x << " " << x_cost << "\n";

        auto insertable = [](Cost cost) { return 0 < cost && cost != invalid_cost; };
        assert(insertable(x_cost));

        // If the subgraph has no neighbors on at least one vertex pair it can be skipped.
        const auto[num_pairs_with_neighbors, num_neighbors_upper_bound] =  m_packing.get_neighbor_count_estimate(x);
        if (num_pairs_with_neighbors < 1 || num_neighbors_upper_bound < 1)
            return {false, false}; // No improvement possible.


        m_packing.add_to_potential(x, x_cost);
        removed_subgraphs.emplace_back(x, x_cost);

#ifndef NDEBUG
        finder->for_all_conversionless_edits(x, [&](auto uv) {
            assert(!m_packing.is_depleted(uv));
            return false;
        });

        VertexPairMap<Cost> potential_copy(m_graph.size());
        for (auto uv : m_graph.vertexPairs())
            potential_copy[uv] = m_packing.potential(uv);
#endif

        // The candidates are partioned by their vertex pairs.
        // For a i (0..pairs.size()-1) the candidates in the range candidates[border[i]..border[i+1]-1] do not contain pairs[0..i-1] and do contain pairs[i].
        auto[pairs, candidates, border] = m_packing.get_closed_neighbors(Subgraph(x));
        // Note: May change to open neighborhood. This is a change of the plateau search, in the sense that the bound is guaranteed to change if an alternative is available.

        // Note: This is true, because of get_closed_neighbors. See the earlier notice for correct solution.
        assert(std::any_of(candidates.begin(), candidates.end(), [&](const auto &c) { return x == c; }));

        // The information about the best solution.
        // Remark: the subgraph x will always be a candidate.
        Cost max_cost = x_cost;
        std::vector<std::vector<std::pair<size_t, Cost>>> max_sets;


        for (size_t pair_i = 0; pair_i < pairs.size(); ++pair_i) {
            for (size_t a_i = border[pair_i]; a_i < border[pair_i + 1]; ++a_i) {
                const Subgraph &a = candidates[a_i];

                auto a_cost = m_packing.calculate_min_cost(a);
                assert(insertable(a_cost));

                m_packing.subtract_from_potential(a, a_cost);

                if (a_cost > max_cost) {
                    max_cost = a_cost;
                    max_sets = {{{a_i, a_cost}}};
                } else if (a_cost == max_cost) {
                    max_sets.push_back({{a_i, a_cost}});
                }

                // Note: The unweighted packing allows the following optimization.
                // `pair_j` can be initialized with `pair_i + 1`. This can be done because the order of (a, b) does not
                // matter and pairs[pair_i] is guaranteed to be blocked by a.
                // This is not the case for weighted packing. pairs[pair_i] is only blocked if it is the pair with the
                // smallest cost. Additionally the order in which the subgraphs are added can change their cost.
                for (size_t pair_j = 0; pair_j < pairs.size(); ++pair_j) {
                    if (m_packing.is_depleted(pairs[pair_j]))
                        continue;
                    for (size_t b_i = border[pair_j]; b_i < border[pair_j + 1]; ++b_i) {
                        if (b_i == a_i) continue;
                        const Subgraph &b = candidates[b_i];

                        auto b_cost = m_packing.calculate_min_cost(b);

                        if (insertable(b_cost)) {
                            // b can be inserted.

                            // Additional linear scan for a third or more subgraphs.
                            // This ensures that each set in max_sets leads to a maximal packing.
                            m_packing.subtract_from_potential(b, b_cost);

                            Cost additional_cost = a_cost + b_cost;
                            std::vector<std::pair<size_t, Cost>> additional_subgraphs;
                            for (size_t pair_k = 0; pair_k < pairs.size(); ++pair_k) {
                                if (m_packing.is_depleted(pairs[pair_k])) // Note: At least two pairs are depleted.
                                    continue;
                                for (size_t c_i = border[pair_k]; c_i < border[pair_k + 1]; ++c_i) {
                                    if (c_i == a_i || c_i == b_i) continue;
                                    const Subgraph &c = candidates[c_i];

                                    auto c_cost = m_packing.calculate_min_cost(c);

                                    if (insertable(c_cost)) {
                                        additional_cost += c_cost;
                                        additional_subgraphs.emplace_back(c_i, c_cost);
                                        m_packing.subtract_from_potential(c, c_cost);
                                    }
                                }
                            }

                            // Remove b and all c's from packing.
                            m_packing.add_to_potential(b, b_cost);
                            for (auto[c_i, c_cost] : additional_subgraphs) {
                                m_packing.add_to_potential(candidates[c_i], c_cost);
                            }
                            // Now only a is in the packing.


                            Cost abc_cost = a_cost + b_cost + additional_cost;
                            additional_subgraphs.emplace_back(a_i, a_cost);
                            additional_subgraphs.emplace_back(b_i, b_cost);

                            if (abc_cost > max_cost) {
                                max_cost = abc_cost;
                                max_sets.clear();
                                max_sets.push_back(std::move(additional_subgraphs));
                            } else if (abc_cost == max_cost) {
                                max_sets.push_back(std::move(additional_subgraphs));
                            }
                        }
                    }
                }

                m_packing.add_to_potential(a, a_cost);
                // No candidates remain in the packing.
            }
        }

#ifndef NDEBUG
        for (auto uv : m_graph.vertexPairs())
            assert(potential_copy[uv] == m_packing.potential(uv));
#endif

        assert(!max_sets.empty());

        std::uniform_int_distribution<size_t> sample(0, max_sets.size() - 1);
        size_t j = sample(m_gen);

        // Must be done here, because candidates will be moved into inserted_subgraphs.
        bool changed = max_cost > x_cost ||
            !(max_sets[j].size() == 1 && candidates[max_sets[j].front().first] == x);

        for (auto[y_i, y_cost] : max_sets[j]) {
            auto &y = candidates[y_i];
            m_packing.subtract_from_potential(y, y_cost);
            inserted_subgraphs.emplace_back(std::move(y), y_cost);
        }

        assert(m_packing.is_valid());
        return {changed, max_cost > x_cost};
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKINGLOCALSEARCH_H
