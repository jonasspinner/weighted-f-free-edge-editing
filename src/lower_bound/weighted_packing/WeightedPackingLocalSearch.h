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

        [[nodiscard]] inline auto &insertion_cost(const Subgraph &subgraph) {
            return m_subgraphs_in_packing[subgraph];
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
            other.reset();

            for (const auto &[subgraph, cost] : m_subgraphs_in_packing) {
                other.subtract_from_potential(subgraph, cost);
            }

            assert(other.cost() == packing.cost());
            for (auto uv : other.depleted_graph().vertex_pairs()) {
                if (other.potential(uv) != packing.potential(uv)) std::cerr << uv << std::endl;
                assert(other.potential(uv) == packing.potential(uv));
                if (other.is_depleted(uv) != packing.is_depleted(uv)) std::cerr << uv << std::endl;
                assert(other.is_depleted(uv) == packing.is_depleted(uv));
            }
        }

#endif

        void reset() {
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

    const EditState *m_edit_state;
    const SubgraphStats<SetOfForbiddenSubgraphs> *m_subgraph_stats;

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
    /**
     * A local search based lower bound algorithm based on a weighted packing.
     *
     * The instance holds a stack of states. Each state keeps track of the forbidden subgraphs which are considered to
     * be in the packing at a certain recursion depth. This is done, as there is no straightforward way to undo the
     * changes (other than maybe by storing a delta between the recursion levels).
     *
     * The weighted packing keeps some kind of residual called potential. Initially the potential is the editing cost
     * for each vertex pair. It is decremented when a incident subgraph is inserted into the packing with a given cost.
     * From the potential we can additional subgraphs which can be inserted into the packing. As the weighted packing
     * itself is quite big with O(n^2) space, there is only one. Most of the times it is consistent with current state,
     * i.e. the state at the deepest recursion level. While ascending from a recursion level, the packing is updated to
     * match the appropriate state.
     *
     * @param instance
     * @param marked
     * @param subgraph_stats
     * @param seed
     */
    WeightedPackingLocalSearch(const EditState *edit_state,
                               const SubgraphStats<SetOfForbiddenSubgraphs> *subgraph_stats,
                               std::size_t seed = 0) :
            m_edit_state(edit_state), m_subgraph_stats(subgraph_stats),
            m_packing(m_edit_state, m_subgraph_stats), m_gen(seed) {
        m_states.emplace_back();
    }

    void push_state(Cost /*k*/) override {
#ifndef NDEBUG
        if (verbosity > 0)
            std::cout << "push_state(...) " << m_states.size() << "\n";
#endif
        if (m_states.empty())
            throw std::runtime_error("There is no state left on the stack. pop_stack has been called to many times.");

        m_states.emplace_back(m_states.back());
    }

    void pop_state() override {
#ifndef NDEBUG
        if (verbosity > 0)
            std::cout << "pop_state() " << m_states.size() << "\n";
#endif
        if (m_states.empty())
            throw std::runtime_error("No state left to pop from the stack.");

        m_states.pop_back();

        // Restore packing from state.
        auto &state = current_state();
        m_packing.reset();
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
     * The vertex pair uv has been marked, but no state has been pushed and no edit has been made.
     *
     * Restore the potential for the vertex pair. An marked vertex pair is ignored while calculating minimum editing
     * costs for forbidden subgraphs. When the potential of the vertex pair is restored to its original editing cost, it
     * may no longer be depleted, i.e. the vertex pair is no longer forbidden when listing subgraphs.
     *
     * Note, that the packing may be no longer maximal. The states will be made maximal in before_edit and after_edit.
     *
     * @param uv
     */
    void after_mark(VertexPair uv) override {
#ifndef NDEBUG
        assert(m_edit_state->is_marked(uv));

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
     * The vertex pair uv has been marked and a state has been pushed, but uv has not been edited yet. The parent state
     * will be unaffected by the coming edit, but uv being marked may have allow more forbidden subgraphs to be inserted
     * into the parent state. This method makes the parent state maximal again.
     *
     * The vertex pair uv is also marked for the current state, but changing the current state does not bring any
     * improvements. The forbidden subgraphs at uv will be either destroyed by the coming edit of uv in the current
     * state, or be converted into another forbidden subgraph. The method `after_edit` is responsible for the current
     * state.
     *
     * The attribute `m_packing` is modified while calculating the changes for the parent state. These change are
     * applied to the parent state, but reverted for `m_packing`. After returning from this method, the packing is the
     * same as it started and consistent with the current state.
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
        if (m_edit_state->cost(uv) == 0) {
            return;
        }

        auto &parent = parent_state();

        auto [solvable, inserted] = insert_incident_subgraphs_into_packing(uv, m_finder, m_edit_state->graph(), m_packing);

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
            parent.insertion_cost(subgraph) += cost;
        }

#ifndef NDEBUG
        current_state().packing_is_same_as_initialized_by_state(m_packing);
#endif
    }


    /**
     * Insert forbidden subgraphs with the vertex pair uv into the packing. No state is modified. If a forbidden
     * subgraph is found, which is fully marked, i.e. cannot be destroyed, the output reflects that the
     * current problem is not solvable any more.
     *
     * @param uv
     * @param finder
     * @param graph
     * @param packing
     *      Will be modified.
     * @return (solvable, inserted_subgraphs)
     */
    [[maybe_unused]] static std::pair<bool, std::vector<std::pair<Cost, Subgraph>>>
    insert_incident_subgraphs_into_packing(VertexPair uv, Finder &finder, const Graph &graph,
            WeightedPacking<SetOfForbiddenSubgraphs> &packing) {
        constexpr bool use_greedy_local_update = false;
        if constexpr (use_greedy_local_update) {
            return insert_incident_subgraphs_into_packing_greedy(uv, finder, graph, packing);
        } else {
            return insert_incident_subgraphs_into_packing_linear(uv, finder, graph, packing);
        }
    }


    [[maybe_unused]] static std::pair<bool, std::vector<std::pair<Cost, Subgraph>>>
    insert_incident_subgraphs_into_packing_linear(VertexPair uv, Finder &finder, const Graph &graph,
            WeightedPacking<SetOfForbiddenSubgraphs> &packing) {
        std::vector<Subgraph> incident_subgraphs;
        std::vector<std::pair<Cost, Subgraph>> inserted;

        finder.find_near(uv, graph, packing.depleted_graph(),
                [&](const Subgraph &subgraph) noexcept {
            assert(subgraph.contains(uv));
            incident_subgraphs.push_back(std::move(subgraph));
            return subgraph_iterators::IterationControl::Continue;
        });

        bool unsolvable = false;
        for (auto &&subgraph : incident_subgraphs) {
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
    insert_incident_subgraphs_into_packing_greedy(VertexPair uv, Finder &finder, const Graph &graph,
            WeightedPacking<SetOfForbiddenSubgraphs> &packing) {
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
     * The vertex pair uv has been marked, a new state has been pushed and uv has been edited.
     *
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
        assert(m_edit_state->is_marked(uv));
        if (verbosity > 0) {
            std::cout << "after_edit(" << uv << ")\n";
            std::cout << "solvable = " << state.is_solvable() << "\n";
        }

        assert(m_packing.is_valid());

        state.packing_is_same_as_initialized_by_state(m_packing);

        assert(m_packing.potential(uv) == m_edit_state->cost(uv));
#endif

        // Remove subgraphs which are destroyed.
        // Note: This can "free" subgraphs which now can be inserted. These are at pairs which were previously depleted.
        std::vector<VertexPair> pairs = remove_incident_subgraphs(uv, m_edit_state->marked_map(), state.subgraphs(), m_packing);
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

    /**
     * Forbidden subgraph with the vertex pair uv are removed, both from map and the packing.
     *
     * @param uv
     * @param marked
     * @param map
     * @param packing
     * @return The vertex pairs which are no longer depleted.
     */
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

    /**
     * Insert subgraphs into the map and the packing.
     *
     * @param subgraphs
     * @param packing
     * @param map
     * @return false if a subgraph was unsolvable, true otherwise.
     */
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

    /**
     * The state and packing are reset and initialized.
     *
     * Precondition:
     *      There is a single state on the stack. This should be the case if there was a balanced number of `push_state`
     *      and `pop_state` calls.
     * @param k
     */
    void initialize(Cost k) override {
        if (verbosity > 0) std::cout << "initialize\n";

        if (m_states.size() != 1)
            throw std::runtime_error("Illegal state. There must be an equal amount of push and pops from the state stack.");

        auto &state = current_state();
        greedy_initialize(state, m_packing, k, m_finder, m_edit_state->graph());

#ifndef NDEBUG
        state.packing_is_same_as_initialized_by_state(m_packing);
#endif
    }

    Cost calculate_lower_bound(Cost k) override {
        if (verbosity > 0) std::cout << "calculate_lower_bound\n";

        auto &state = current_state();

        if (state.is_solvable()) {

#ifndef NDEBUG
            assert(m_packing.is_maximal());
            assert(m_packing.is_valid());
            state.packing_is_same_as_initialized_by_state(m_packing);
#endif

            if (m_packing.cost() <= k)
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

#ifndef NDEBUG
            assert(m_packing.is_maximal());
            assert(m_packing.is_valid());
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

            for (const auto &[subgraph, cost] : state.subgraphs()) {
                auto[changed, improved] = \
                    find_one_two_improvement(subgraph, cost, m_removed_subgraphs, m_inserted_subgraphs);
                changed_in_round |= changed;
                improved_in_round |= improved;
            }

            // Update state. During a round, the state is stale and updates are only represented in m_packing.
            // Insert before remove, because both could cancel each other.
            for (auto&&[subgraph, cost] : m_inserted_subgraphs) {
                // If the subgraph is not yet in the map, the initial value is 0.
                state.insertion_cost(subgraph) += cost;
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
            assert(m_packing.is_maximal());
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

    void greedy_initialize(State &state, WeightedPacking<SetOfForbiddenSubgraphs> &packing, Cost k, Finder& finder,
                           const Graph& graph) const {
#ifndef NDEBUG
        if (verbosity > 0) std::cout << "greedy_initialize\n";
#endif

        state.reset();
        packing.reset();

        // Note: Could turn this into an attribute and call clear() at beginning. This would prevent reallocating memory.
        std::vector<std::pair<Subgraph, Cost>> subgraph_heap;

        auto comp = [](const auto &lhs, const auto &rhs) { return lhs.second < rhs.second; };

        auto exit_state = finder.find(graph, packing.depleted_graph(), [&](const Subgraph &subgraph) {
            Cost initial_min_cost = packing.calculate_min_cost(subgraph);

            subgraph_heap.emplace_back(subgraph, initial_min_cost);

            return subgraph_iterators::break_if(initial_min_cost > k);
        });

        if (exit_state == subgraph_iterators::IterationExit::Break) {
            state.set_unsolvable();
            return;
        }

        std::make_heap(subgraph_heap.begin(), subgraph_heap.end(), comp);

        while (!subgraph_heap.empty() && packing.cost() <= k) {
            std::pop_heap(subgraph_heap.begin(), subgraph_heap.end(), comp);
            auto &[subgraph, cost] = subgraph_heap.back();

            Cost current_min_cost = packing.calculate_min_cost(subgraph);

            if (current_min_cost == 0) {
                // The subgraph cannot be inserted.
                subgraph_heap.pop_back();
            } else if (current_min_cost == cost) {
                // The editing cost is maximal from all remaining subgraphs in the queue.
                // Update packing and state.
                packing.subtract_from_potential(subgraph, cost);
                state.insertion_cost(subgraph) = cost;

                subgraph_heap.pop_back();
            } else {
                // The editing cost is no longer up to date. The subgraph may be inserted in the future.
                subgraph_heap.back().second = current_min_cost;
                std::push_heap(subgraph_heap.begin(), subgraph_heap.end(), comp);
            }
        }

        if (m_packing.cost() <= k) {
            assert(packing.is_maximal());
        } else {
            state.set_unsolvable();
        }
    }

    /**
     *  Find (1, 2) improvements.
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
    std::tuple<bool, bool> find_one_two_improvement(const Subgraph &x, Cost x_cost,
                                                     std::vector<std::pair<Subgraph, Cost>> &removed_subgraphs,
                                                     std::vector<std::pair<Subgraph, Cost>> &inserted_subgraphs) {
#ifndef NDEBUG
        if (verbosity > 0) std::cout << "find_one_two_improvement " << x << " " << x_cost << "\n";
#endif

        assert(x_cost >= 1);
        m_packing.add_to_potential(x, 1);
        auto[pairs, candidates, border] = m_packing.get_open_neighbors(x, 1);

#ifndef NDEBUG
        if (verbosity > 0) std::cout << "candidates.size() = " << candidates.size() << "\n";
#endif

        if (candidates.empty()) {
            m_packing.subtract_from_potential(x, 1);
            return {false, false};
        }


        auto insertable = [&](const Subgraph &subgraph) {
            return m_packing.calculate_min_cost(subgraph) > 0;
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
            assert(cost > 0);
            assert(cost != invalid_cost);
#ifndef NDEBUG
            if (verbosity > 0)
                std::cout << "inserted " << subgraph << " " << cost << "\n";
#endif
            m_packing.subtract_from_potential(subgraph, cost);
            inserted_subgraphs.emplace_back(subgraph, cost);
        };

        auto try_insert = [&](const auto &subgraph) {
            auto cost = m_packing.calculate_min_cost(subgraph);
            assert(cost != invalid_cost);
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


        if (x_cost > 1) {
            m_packing.add_to_potential(x, x_cost - 1);  // Remaining cost, as x was already removed with cost 1.
        }
        removed_subgraphs.emplace_back(x, x_cost);

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


    /**
     * Select a subgraph from the candidates and the additional subgraph. This is similar to a Maximum Weight
     * Independent Set heuristic.
     *
     * @param candidates
     * @param subgraph
     * @return
     */
    Subgraph find_plateau_candidate(const std::vector<Subgraph>& candidates, const Subgraph &subgraph) {
        constexpr auto plateau_search_use_min_degree = true;
        if (plateau_search_use_min_degree) {
            return find_plateau_candidate_min_degree(candidates, subgraph, *m_subgraph_stats);
        } else {
            return find_plateau_candidate_random(candidates, m_gen);
        }
    }

    static Subgraph find_plateau_candidate_random(const std::vector<Subgraph>& candidates,
            std::mt19937_64 &gen) {
        std::uniform_int_distribution<size_t> dist(0, candidates.size() - 1);
        return candidates[dist(gen)];
    }

    static Subgraph find_plateau_candidate_min_degree(const std::vector<Subgraph>& candidates,
            const Subgraph & subgraph, const SubgraphStats<SetOfForbiddenSubgraphs> &subgraph_stats) {
        auto estimate_degree = [&](auto &s) {
            std::size_t degree_ub = 0;
            for (auto uv : s.non_converting_edits()) {
                degree_ub += subgraph_stats.subgraph_count(uv);
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
