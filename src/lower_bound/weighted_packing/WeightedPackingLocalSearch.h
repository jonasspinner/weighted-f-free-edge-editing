#ifndef WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKINGLOCALSEARCH_H
#define WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKINGLOCALSEARCH_H

#include "WeightedPacking.h"

#include "robin_hood.h"
#include <unordered_set>


namespace lower_bound {

template <Options::FSG SetOfForbiddenSubgraphs>
class WeightedPackingLocalSearch final : public LowerBoundI {
    using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
    using Finder = typename Subgraph::Finder;

    class State {
    public:
        using Map = robin_hood::unordered_map<Subgraph, Cost>;
    private:
        Map m_subgraphs_in_packing;
        bool m_solvable = true;
    public:
        [[nodiscard]] inline auto &subgraphs() {
            return m_subgraphs_in_packing;
        }

        [[nodiscard]] inline const auto &subgraphs() const {
            return m_subgraphs_in_packing;
        }

        void set_unsolvable() {
            m_solvable = false;
        }

        [[nodiscard]] bool is_solvable() const {
            return m_solvable;
        }

#ifndef NDEBUG

        void packing_is_same_as_initialized_by_state(const WeightedPacking<SetOfForbiddenSubgraphs> &packing) {
            WeightedPacking<SetOfForbiddenSubgraphs> other(packing);
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
    const SubgraphStats<SetOfForbiddenSubgraphs> &m_subgraph_stats;

    WeightedPacking<SetOfForbiddenSubgraphs> m_packing;
    std::vector<State> m_states;

    std::mt19937_64 m_gen;

#ifndef NDEBUG
    int verbosity = 0;
#else
    static constexpr int verbosity = 0;
#endif

    Finder m_finder;

    std::vector<std::pair<Subgraph, Cost>> m_removed_subgraphs, m_inserted_subgraphs;
public:
    WeightedPackingLocalSearch(const Instance &instance, const VertexPairMap<bool> &marked,
                               const SubgraphStats<SetOfForbiddenSubgraphs> &subgraph_stats,
                               std::size_t seed = 0) :
            m_graph(instance.graph), m_costs(instance.costs), m_marked(marked), m_subgraph_stats(subgraph_stats),
            m_packing(instance, marked, subgraph_stats), m_gen(seed) {
        m_states.emplace_back();
    }

    void push_state(Cost /*k*/) override {
#ifndef NDEBUG
        if (verbosity > 0)
            std::cout << "push_state(...) " << m_states.size() << "\n";
#endif
        assert(!m_states.empty());
        m_states.emplace_back(m_states.back());
    }

    void pop_state() override {
#ifndef NDEBUG
        if (verbosity > 0)
            std::cout << "pop_state() " << m_states.size() << "\n";
#endif
        assert(!m_states.empty());
        m_states.pop_back();

        // Restore packing from state.
        auto &state = current_state();
        m_packing.clear();
        for (const auto &[subgraph, cost] : state.subgraphs()) {
            m_packing.subtract_from_potential(subgraph, cost);
        }
    }

    inline State &current_state() {
        assert(!m_states.empty());
        return m_states.back();
    }

    inline State &parent_state() {
        assert(m_states.size() > 1);
        return m_states[m_states.size() - 2];
    }

    /**
     * Restore the potential for the vertex pair. When a vertex pair is marked, it is ignored when the editing costs for
     * subgraphs are calculated.
     *
     * Note, that the packing may be no longer maximal. The states will be made maximal in before_edit and after_edit.
     *
     * @param uv
     */
    void after_mark(VertexPair uv) override {
#ifndef NDEBUG
        assert(m_marked[uv]);

        if (verbosity > 0)
            std::cout << "after_mark(" << uv << ")\n";

        auto &state = current_state();

        assert(m_packing.is_valid());
        if (state.is_solvable()) {
            assert(m_packing.is_maximal());
        }
#endif

        m_packing.restore_potential(uv);

#ifndef NDEBUG
        state.packing_is_same_as_initialized_by_state(m_packing);
#endif
    }

    /**
     * Make parent state maximal.
     *
     * m_packing is used to calculate the changes for the parent state. These change are applied to the parent state,
     * but undone for m_packing. After the method, the packing is the same as it started as.
     *
     * Precondition:
     *      current and parent state are identical
     *
     * Postcondititon:
     *      if parent state is solvable, then it is maximal
     *
     * Invariant:
     *      m_packing is consistent with the current state
     *
     * @param uv
     */
    void before_edit(VertexPair uv) override {
#ifndef NDEBUG
        current_state().packing_is_same_as_initialized_by_state(m_packing);
#endif
        if (m_costs[uv] == 0) {
            return;
        }

        auto &parent = parent_state();

        auto [solvable, inserted] = insert_incident_subgraphs_into_packing(uv, m_finder, m_graph, m_packing);

        // Packing reflects parent_state() after it has been updated. It should be maximal or unsolvable.
        if (solvable) {
            assert(m_packing.is_maximal());
        } else {
            parent.set_unsolvable();
        }


        for (auto &[cost, subgraph] : inserted) {
            // Undo changes in packing to reflect current_state().
            m_packing.add_to_potential(subgraph, cost);

            // Commit insertions to parent_state();
            parent.subgraphs()[subgraph] += cost;
        }

#ifndef NDEBUG
        current_state().packing_is_same_as_initialized_by_state(m_packing);
#endif
    }


    [[maybe_unused]] static std::pair<bool, std::vector<std::pair<Cost, Subgraph>>>
    insert_incident_subgraphs_into_packing(VertexPair uv, Finder &finder, const Graph &graph, WeightedPacking<SetOfForbiddenSubgraphs> &packing) {
        constexpr bool use_greedy_local_update = false;
        if constexpr (use_greedy_local_update) {
            return insert_incident_subgraphs_into_packing_greedy(uv, finder, graph, packing);
        } else {
            return insert_incident_subgraphs_into_packing_linear(uv, finder, graph, packing);
        }
    }


    [[maybe_unused]] static std::pair<bool, std::vector<std::pair<Cost, Subgraph>>>
    insert_incident_subgraphs_into_packing_linear(VertexPair uv, Finder &finder, const Graph &graph, WeightedPacking<SetOfForbiddenSubgraphs> &packing) {
        std::vector<Subgraph> incident_subgraphs;
        std::vector<std::pair<Cost, Subgraph>> inserted;

        finder.find_near(uv, graph, packing.depleted_graph(),
                [&](const Subgraph &subgraph) noexcept {
            assert(subgraph.contains(uv));
            incident_subgraphs.push_back(std::move(subgraph));
            return false;
        });

        bool unsolvable = false;
        for (auto&& subgraph : incident_subgraphs) {
            auto cost = packing.calculate_min_cost(subgraph);
            if (cost == 0) continue;
            if (cost == invalid_cost) {
                unsolvable = true;
                break;
            }
            packing.subtract_from_potential(subgraph, cost);
            inserted.emplace_back(cost, std::move(subgraph));
        }

        return {(!unsolvable), std::move(inserted)};
    }


    [[maybe_unused]] static std::pair<bool, std::vector<std::pair<Cost, Subgraph>>>
    insert_incident_subgraphs_into_packing_greedy(VertexPair uv, Finder &finder, const Graph &graph, WeightedPacking<SetOfForbiddenSubgraphs> &packing) {
        std::vector<std::pair<Cost, Subgraph>> subgraph_heap;

        bool unsolvable = finder.find_near(uv, graph, packing.depleted_graph(),
                [&](const Subgraph &subgraph) noexcept {
            assert(subgraph.contains(uv));
            auto cost = packing.calculate_min_cost(subgraph);
            if (cost == 0) return false;
            if (cost == invalid_cost) {
                return true;
            }
            subgraph_heap.emplace_back(cost, subgraph);
            return false;
        });

        if (unsolvable) {
            return {false, {}};
        }

        // std::cout << "before_edit subgraph_heap.size() = " << subgraph_heap.size() << "\n";
        std::make_heap(subgraph_heap.begin(), subgraph_heap.end());
        std::vector<std::pair<Cost, Subgraph>> inserted;

        while (!subgraph_heap.empty()) {
            std::pop_heap(subgraph_heap.begin(), subgraph_heap.end());
            const auto &[old_cost, subgraph] = subgraph_heap.back();

            // m_cost_remaining may be updated after cost has been calculated.
            Cost current_cost = packing.calculate_min_cost(subgraph);

            if (current_cost == 0) {
                // The subgraph cannot be inserted.
                subgraph_heap.pop_back();
            } else if (current_cost == old_cost) {
                // Subgraph can be inserted.
                // Update packing to keep min cost calculation consistent. The subgraphs will actually be inserted later
                // and the changes to the packing will be undone.
                packing.subtract_from_potential(subgraph, current_cost);
                inserted.push_back(std::move(subgraph_heap.back()));
                subgraph_heap.pop_back();
            } else {
                // The editing cost is no longer up to date. The subgraph may be inserted in the future.
                subgraph_heap.back().first = current_cost;
                std::push_heap(subgraph_heap.begin(), subgraph_heap.end());
            }
        }

        return {true, std::move(inserted)};
    }


    /**
     * Make current state maximal. Incident subgraphs may have been destroyed and are no longer valid. Replace them by
     * newly created subgraphs.
     *
     * Precondiditions:
     *      m_packing is consistent with the current state
     *      uv is marked and its potential is restored
     *
     * Postconditions:
     *      if current state becomes unsolvable, the state of m_packing is undefined
     *      else
     *      m_packing is consistent with the current state
     *      m_packing is maximal
     *
     * @param uv
     */
    void after_edit(VertexPair uv) override {

        auto &state = current_state();

#ifndef NDEBUG
        assert(m_marked[uv]);
        if (verbosity > 0) {
            std::cout << "after_edit(" << uv << ")\n";
            std::cout << "solvable = " << state.is_solvable() << "\n";
        }

        assert(m_packing.is_valid());

        state.packing_is_same_as_initialized_by_state(m_packing);

        assert(m_packing.potential(uv) == m_costs[uv]);
#endif

        // Remove subgraphs which are destroyed.
        // Note: This can "free" subgraphs which now can be inserted. These are at pairs which were previously depleted.
        std::vector<VertexPair> pairs = remove_incident_subgraphs(uv, m_marked, state.subgraphs(), m_packing);
        if (!m_packing.is_depleted(uv))
            pairs.push_back(uv);


#ifndef NDEBUG
        if (verbosity > 0) {
            if (!pairs.empty())
                std::cout << "### pairs :";
            for (auto xy : pairs)
                std::cout << " " << xy;
            std::cout << std::endl;
        }

        for (auto xy : pairs) {
            assert(!m_packing.is_depleted(xy));
        }
#endif

        auto subgraphs = m_packing.get_incident_subgraphs(pairs);

        bool solvable = insert_subgraphs(subgraphs, m_packing, state.subgraphs());

        if (!solvable) {
            state.set_unsolvable();
            return;
        }

#ifndef NDEBUG
        assert(m_packing.is_maximal());
        assert(m_packing.is_valid());

        state.packing_is_same_as_initialized_by_state(m_packing);
#endif
    }

    static std::vector<VertexPair> remove_incident_subgraphs(VertexPair uv, const VertexPairMap<bool> &marked,
            typename State::Map &map, WeightedPacking<SetOfForbiddenSubgraphs> &packing) {
        std::vector<VertexPair> pairs;

        for (auto it = map.begin(); it != map.end();) {
            const auto &[subgraph, cost] = *it;
            auto edits = subgraph.non_converting_edits();
            if (std::any_of(edits.begin(), edits.end(), [&](auto xy) { return xy == uv; })) {
                for (auto xy : edits) {
                    if (!marked[xy] && packing.is_depleted(xy)) {
                        pairs.push_back(xy);
                    }
                }

                packing.add_to_potential(subgraph, cost);
                it = map.erase(it);
            } else {
                ++it;
            }
        }

        return pairs;
    }

    static bool insert_subgraphs(std::vector<Subgraph> &subgraphs,
            WeightedPacking<SetOfForbiddenSubgraphs> &packing, typename State::Map &map) {
        constexpr bool use_greedy_insertions = false;

        if constexpr (use_greedy_insertions) {
            return insert_subgraphs_into_packing_and_map_greedy(subgraphs, packing, map);
        } else {
            return insert_subgraphs_into_packing_and_map_linear(subgraphs, packing, map);
        }
    }

    [[maybe_unused]] static bool insert_subgraphs_into_packing_and_map_linear(std::vector<Subgraph> &subgraphs,
            WeightedPacking<SetOfForbiddenSubgraphs> &packing, typename State::Map &map) {

        for (auto &subgraph : subgraphs) {
            auto cost = packing.calculate_min_cost(subgraph);
            if (cost == invalid_cost) {
                return false;
            } else if (cost > 0) {
                packing.subtract_from_potential(subgraph, cost);
                map[std::move(subgraph)] += cost;
            }
        }
        return true;
    }

    [[maybe_unused]] static bool insert_subgraphs_into_packing_and_map_greedy(std::vector<Subgraph> &subgraphs,
            WeightedPacking<SetOfForbiddenSubgraphs> &packing, typename State::Map &map) {
        std::vector<std::pair<Cost, Subgraph>> subgraph_heap;

        for (auto &subgraph : subgraphs) {
            auto cost = packing.calculate_min_cost(subgraph);
            if (cost == invalid_cost) {
                return false;
            }
            subgraph_heap.emplace_back(cost, std::move(subgraph));
        }

        // std::cout << "after_edit subgraph_heap.size() = " << subgraph_heap.size() << "\n";
        std::make_heap(subgraph_heap.begin(), subgraph_heap.end());

        while (!subgraph_heap.empty()) {
            std::pop_heap(subgraph_heap.begin(), subgraph_heap.end());
            auto &[old_cost, subgraph] = subgraph_heap.back();

            // m_cost_remaining may be updated after cost has been calculated.
            Cost current_cost = packing.calculate_min_cost(subgraph);

            if (current_cost == 0) {
                // The subgraph cannot be inserted.
                subgraph_heap.pop_back();
            } else if (current_cost == old_cost) {
                // Subgraph can be inserted.
                packing.subtract_from_potential(subgraph, current_cost);
                map[std::move(subgraph)] += current_cost;
                subgraph_heap.pop_back();
            } else {
                // The editing cost is no longer up to date. The subgraph may be inserted in the future.
                subgraph_heap.back().first = current_cost;
                std::push_heap(subgraph_heap.begin(), subgraph_heap.end());
            }
        }

        return true;
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
            // auto before = m_packing.cost();
            if (m_packing.cost() <= k)
                local_search(state, k);
            // std::cout << "# " << k << " " << before << " " << (m_packing.cost() - before) << " " << (k - m_packing.cost()) << "\n";
            return m_packing.cost();
        } else {
            // std::cout << k << "unsolvable\n";
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

        size_t max_iter = std::numeric_limits<std::size_t>::max();
        size_t max_rounds_no_improvement = 5;
        size_t num_rounds_no_improvement = 0;

        if (state.subgraphs().empty())
            return;

        assert(m_packing.is_maximal());

        for (size_t iter = 0; iter < max_iter; ++iter) {
            bool changed_in_round = false, improved_in_round = false;
            m_removed_subgraphs.clear();
            m_inserted_subgraphs.clear();

            // std::cout << iter << " " << state.subgraphs().size() << " " << m_packing.cost() << "\n";
            for (const auto &[subgraph, cost] : state.subgraphs()) {
                auto[changed, improved] = \
                    find_one_two_improvement_2(subgraph, cost, m_removed_subgraphs, m_inserted_subgraphs);
                changed_in_round |= changed;
                improved_in_round |= improved;
            }

            // Update state. During a round, the state is stale and updates are only represented in m_packing.
            // Insert before remove, because both could cancel each other.
            for (auto&&[subgraph, cost] : m_inserted_subgraphs) {
                // If the subgraph is not yet in the map, the initial value is 0.
                state.subgraphs()[std::move(subgraph)] += cost;
            }
            for (auto&&[subgraph, cost] : m_removed_subgraphs) {
                // If the subgraph is fully removed, i.e. its cost is zero, it is removed from the map.
                auto it = state.subgraphs().find(subgraph);
                assert(it->second >= cost);
                it->second -= cost;
                if (it->second == 0) {
                    state.subgraphs().erase(it);
                }
            }

#ifndef NDEBUG
            state.packing_is_same_as_initialized_by_state(m_packing);
            assert(m_packing.is_maximal()); // TODO: fix
#endif

            num_rounds_no_improvement = improved_in_round ? 0 : num_rounds_no_improvement + 1;
            if (iter >= max_iter ||
                !changed_in_round ||
                num_rounds_no_improvement > max_rounds_no_improvement ||
                m_packing.cost() > k) {
//                std::cout << iter + 1 << " ";
//                std::cout << !changed_in_round << " ";
//                std::cout << num_rounds_no_improvement;
//                if (m_packing.cost() > k)
//                    std::cout << " bound";
//                std::cout << "\n";
                break;
            }
        }
    }

    void greedy_initialize(State &state, Cost k) {
#ifndef NDEBUG
        if (verbosity > 0)
            std::cout << "greedy_initialize\n";
#endif
        std::vector<std::pair<Subgraph, Cost>> subgraph_heap; // Note: Could turn this into an attribute and call clear() at beginning. This would prevent reallocating memory.

        auto comp = [](const auto &lhs, const auto &rhs) { return lhs.second < rhs.second; };

        bool unsolvable = m_finder.find(m_graph, m_packing.depleted_graph(), [&](const Subgraph &subgraph) {
            Cost initial_min_cost = m_packing.calculate_min_cost(subgraph);

            subgraph_heap.emplace_back(subgraph, initial_min_cost);

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
        m_finder->for_all_conversionless_edits(x, [&](auto uv) {
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
        // Note: Improvement: only list candidates incident to vertex pairs which became non-depleted by removing x.

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

                for (size_t pair_j = pair_i; pair_j < pairs.size(); ++pair_j) {
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

                            Cost additional_cost = 0;
                            std::vector<std::pair<size_t, Cost>> additional_subgraphs;
                            for (size_t pair_k = 0; pair_k < pairs.size(); ++pair_k) {
                                if (m_packing.is_depleted(pairs[pair_k]))
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

#ifndef NDEBUG
        Cost sum = 0;
        for (auto[y_i, y_cost] : max_sets[j]) {
            sum += y_cost;
        }
        assert(max_cost == sum);
#endif

        // Must be done here, because candidates will be moved into inserted_subgraphs.
        bool changed = max_cost > x_cost ||
            !(max_sets[j].size() == 1 && candidates[max_sets[j].front().first] == x);

        for (auto[y_i, y_cost] : max_sets[j]) {
            auto &y = candidates[y_i];
            m_packing.subtract_from_potential(y, y_cost);
            inserted_subgraphs.emplace_back(std::move(y), y_cost);
        }

        assert(m_packing.is_valid());
        return std::make_pair(changed, max_cost > x_cost);
    }

    /**
     * Alternative version.
     *
     *  + Remove x with cost 1.
     *  + List candidates which are incident to vertex pairs which became non-depleted by removing x. Do not list x as
     *    candidate.
     *  + Iterate over candidates a.
     *  + If another candidate b also fits into a packing, then (a, b) is an improvement.
     *  + Fully remove x.
     *  + Insert (a, b).
     *  + Make packing maximal by inserting additional candidates c_1, ...
     *  + Insert x again.
     *
     *  If no subgraph b has been found, every candidate is a plateau-search candidate.
     *  + Use an heuristic to choose one of the candidates. For example:
     *    + A subgraph which is already in the packing.
     *    + A subgraph which does deplete the fewest vertex pairs.
     *
     * @param x
     * @param x_cost
     * @param removed_subgraphs
     * @param inserted_subgraphs
     * @return
     */
    std::tuple<bool, bool> find_one_two_improvement_2(const Subgraph &x, Cost x_cost,
                                                     std::vector<std::pair<Subgraph, Cost>> &removed_subgraphs,
                                                     std::vector<std::pair<Subgraph, Cost>> &inserted_subgraphs) {
#ifndef NDEBUG
        if (verbosity > 0)
            std::cout << "find_one_two_improvement_2 " << x << " " << x_cost << "\n";
#endif

        assert(x_cost >= 1);
        m_packing.add_to_potential(x, 1);
        auto[pairs, candidates, border] = m_packing.get_open_neighbors(x, 1);

#ifndef NDEBUG
        if (verbosity > 0)
            std::cout << "candidates.size() = " << candidates.size() << "\n";
#endif
        if (candidates.empty()) {
            m_packing.subtract_from_potential(x, 1);
            return {false, false};
        }


        auto insertable = [&](const Subgraph &subgraph) {
            return m_packing.calculate_min_cost(subgraph) > 0;
//            return !m_finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
//                return !m_marked[uv] && m_packing.is_depleted(uv);
//            });
        };


        constexpr auto invalid_index = std::numeric_limits<std::size_t>::max();
        std::size_t a_found = invalid_index;
        std::size_t b_found = invalid_index;

        bool improvement_found = false;

        for (size_t pair_i = 0; pair_i < pairs.size() && !improvement_found; ++pair_i) {
            for (size_t a_i = border[pair_i]; a_i < border[pair_i + 1] && !improvement_found; ++a_i) {
                const Subgraph &a = candidates[a_i];

#ifndef NDEBUG
                assert(insertable(a));
                auto a_cost = m_packing.calculate_min_cost(a);
                assert(a_cost == 1);
#endif

                m_packing.subtract_from_potential(a, 1);

                for (size_t pair_j = pair_i + 1; pair_j < pairs.size() && !improvement_found; ++pair_j) {
                    if (m_packing.is_depleted(pairs[pair_j]))
                        continue;
                    for (size_t b_i = border[pair_j]; b_i < border[pair_j + 1] && !improvement_found; ++b_i) {
                        const Subgraph &b = candidates[b_i];

                        if (insertable(b)) {
                            a_found = a_i;
                            b_found = b_i;
                            improvement_found = true;
                        }
                    }
                }

                m_packing.add_to_potential(a, 1);
            }
        }




        auto insert = [&](const auto &subgraph) {
            auto cost = m_packing.calculate_min_cost(subgraph);
#ifndef NDEBUG
            if (verbosity > 0)
                std::cout << "inserted " << subgraph << " " << cost << "\n";
#endif
            m_packing.subtract_from_potential(subgraph, cost);
            inserted_subgraphs.emplace_back(subgraph, cost);
        };

        auto try_insert = [&](const auto &subgraph) {
            auto cost = m_packing.calculate_min_cost(subgraph);
            if (cost > 0) {
#ifndef NDEBUG
                if (verbosity > 0)
                    std::cout << "inserted " << subgraph << " " << cost << "\n";
#endif
                m_packing.subtract_from_potential(subgraph, cost);
                inserted_subgraphs.emplace_back(subgraph, cost);
                return true;
            }
            return false;
        };


        if (x_cost == 1) {
            removed_subgraphs.emplace_back(x, 1);
        } else {
            assert(x_cost >= 1);
            m_packing.add_to_potential(x, x_cost - 1);
            removed_subgraphs.emplace_back(x, x_cost);
        }

        if (improvement_found) {
#ifndef NDEBUG
            if (verbosity > 0)
                std::cout << "(1,2) improvement\n";
#endif
            assert(a_found != invalid_index);
            assert(b_found != invalid_index);

            // insert(candidates[a_found]);
            m_packing.subtract_from_potential(candidates[a_found], 1);
            inserted_subgraphs.emplace_back(candidates[a_found], 1);

            insert(candidates[b_found]);
            try_insert(candidates[a_found]);
        } else {
#ifndef NDEBUG
            if (verbosity > 0)
                std::cout << "(1,1) plateau\n";
#endif
            auto candidate = find_plateau_candidate(candidates, x);
            insert(candidate);

            if (candidate == x)
                return {false, false};
        }

        for (auto &c : candidates) {
            try_insert(c);
        }

        try_insert(x);

        assert(m_packing.is_maximal());

        return {true, improvement_found};
    }

    Subgraph find_plateau_candidate(const std::vector<Subgraph>& candidates, const Subgraph &subgraph) {
        constexpr auto plateau_search_use_min_degree = true;
        if (plateau_search_use_min_degree) {
            return find_plateau_candidate_min_degree(candidates, subgraph);
        } else {
            return find_plateau_candidate_random(candidates);
        }
    }

    Subgraph find_plateau_candidate_random(const std::vector<Subgraph>& candidates) {
        std::uniform_int_distribution<size_t> dist(0, candidates.size() - 1);
        return candidates[dist(m_gen)];
    }

    Subgraph find_plateau_candidate_min_degree(const std::vector<Subgraph>& candidates, const Subgraph & subgraph) {
        auto estimate_degree = [&](auto &s) {
            std::size_t degree_ub = 0;
            for (auto uv : s.non_converting_edits()) {
                degree_ub += m_subgraph_stats.subgraphCount(uv);
            }
            return degree_ub;
        };

        std::size_t min_degree = estimate_degree(subgraph);
        Subgraph min_degree_subgraph = subgraph;


        for (const auto& c : candidates) {
            auto degree_ub = estimate_degree(subgraph);
            if (degree_ub < min_degree) {
                min_degree = degree_ub;
                min_degree_subgraph = c;
            }
        }
        return min_degree_subgraph;
    }
};

}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKINGLOCALSEARCH_H
