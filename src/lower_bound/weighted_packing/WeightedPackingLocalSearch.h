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

        void set_unsolvable() {
            m_solvable = false;
        }

        [[nodiscard]] bool is_solvable() const {
            return m_solvable;
        }
    };

    const Graph &m_graph;
    const VertexPairMap<Cost> &m_costs;
    const VertexPairMap<bool> &m_marked;

    WeightedPacking m_packing;
    std::vector<std::unique_ptr<State>> m_states;

    std::mt19937_64 m_gen;

public:
    WeightedPackingLocalSearch(const Instance &instance, const VertexPairMap<bool> &marked,
                               const SubgraphStats &subgraph_stats, std::shared_ptr<FinderI> finder_ref) :
            LowerBoundI(std::move(finder_ref)), m_graph(instance.graph), m_costs(instance.costs), m_marked(marked),
            m_packing(instance, marked, subgraph_stats) {
        m_gen.seed(0);
        m_states.emplace_back();
    }

    State &current_state() {
        assert(!m_states.empty());
        return *m_states.back();
    }

    void initialize(Cost k) override {
        auto &state = current_state();
        greedy_initialize(state, k);
    }

    Cost calculate_lower_bound(Cost k) override {
        auto &state = current_state();
        state.subgraphs().clear();
        m_packing.clear();
        greedy_initialize(state, k);
        if (state.is_solvable()) {
            local_search(state);
            return m_packing.cost();
        } else {
            return invalid_cost;
        }
    }

    void local_search(State &state) {
        size_t max_iter = 5;
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
            if (!changed_in_round || num_rounds_no_improvement > max_rounds_no_improvement) {
                break;
            }
        }
    }

    void greedy_initialize(State &state, Cost k) {
        std::vector<std::pair<Subgraph, Cost>> subgraph_heap; // Note: Could turn this into an attribute and call clear() at beginning. This would prevent reallocating memory.
        Cost lower_bound = 0;
        Cost max_min_cost = std::numeric_limits<Cost>::min();

        bool unsolvable = finder->find_with_duplicates(m_graph, [&](Subgraph &&subgraph) {
            Cost initial_min_cost = m_packing.calculate_min_cost(subgraph);

            subgraph_heap.emplace_back(std::move(subgraph), initial_min_cost);

            max_min_cost = std::max(max_min_cost, initial_min_cost);
            return max_min_cost > k;
        });

        if (unsolvable) {
            state.set_unsolvable();
        }

        std::make_heap(subgraph_heap.begin(), subgraph_heap.end());

        while (!subgraph_heap.empty() && lower_bound < k) {
            std::pop_heap(subgraph_heap.begin(), subgraph_heap.end());
            const auto &[subgraph, cost] = subgraph_heap.back();

            Cost current_min_cost = m_packing.calculate_min_cost(subgraph);

            if (current_min_cost == 0) {
                // The subgraph cannot be inserted.
                subgraph_heap.pop_back();
            } else if (current_min_cost == cost) {
                // The editing cost is maximal from all remaining subgraphs in the queue. Increase the lower bound
                // and update the remaining cost matrix.
                lower_bound += current_min_cost;

                m_packing.subtract_from_potential(subgraph, cost);
                state.subgraphs().emplace(std::move(subgraph_heap.back()));

                subgraph_heap.pop_back();
            } else {
                // The editing cost is no longer up to date. The subgraph may be inserted in the future.
                subgraph_heap.back().second = current_min_cost;
                std::push_heap(subgraph_heap.begin(), subgraph_heap.end());
            }
        }

        assert(m_packing.is_maximal());
    }


    std::tuple<bool, bool> find_one_two_improvement(const Subgraph &x, Cost x_cost,
                                                    std::vector<std::pair<Subgraph, Cost>> &removed_subgraphs,
                                                    std::vector<std::pair<Subgraph, Cost>> &inserted_subgraphs) {

        auto insertable = [](Cost cost) { return 0 < cost && cost != invalid_cost; };
        assert(insertable(x_cost));

        // If the subgraph has no neighbors on at least one vertex pair it can be skipped.
        const auto[num_pairs_with_neighbors, num_neighbors_upper_bound] =  m_packing.get_neighbor_count_estimate(x);
        if (num_pairs_with_neighbors < 1 || num_neighbors_upper_bound < 1)
            return {false, false}; // No improvement possible.


        m_packing.add_to_potential(x, x_cost);
        removed_subgraphs.emplace_back(x, x_cost);

        // The candidates are partioned by their vertex pairs.
        // For a i (0..pairs.size()-1) the candidates in the range candidates[border[i]..border[i+1]-1] do not contain pairs[0..i-1] and do contain pairs[i].
        const auto pairs = m_packing.get_neighbor_pairs(x);
        auto[candidates, border] = m_packing.get_neighbors(pairs);

        // The information about the best solution.
        // Remark: the subgraph x will always be a candidate.
        Cost max_cost = x_cost;
        std::vector<std::vector<std::pair<size_t, Cost>>> max_i;


        for (size_t pair_i = 0; pair_i < pairs.size(); ++pair_i) {
            for (size_t a_i = border[pair_i]; a_i < border[pair_i + 1]; ++a_i) {
                const Subgraph &a = candidates[a_i];

                auto a_cost = m_packing.calculate_min_cost(a);
                m_packing.subtract_from_potential(a, a_cost);

                assert(insertable(a_cost));

                if (a_cost > max_cost) {
                    max_cost = a_cost;
                    max_i = {{{a_i, a_cost}}};
                } else if (a_cost == max_cost) {
                    max_i.push_back({{a_i, a_cost}});
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
                        const Subgraph &b = candidates[b_i];

                        auto b_cost = m_packing.calculate_min_cost(b);
                        if (insertable(b_cost)) {
                            // b can be inserted.

                            // Additional linear scan for a third or more subgraphs.
                            // This ensures that each set in max_i leads to a maximal packing.
                            m_packing.subtract_from_potential(b, b_cost);

                            Cost additional_cost = a_cost + b_cost;
                            std::vector<std::pair<size_t, Cost>> additional_subgraphs;
                            for (size_t pair_k = 0; pair_k < pairs.size(); ++pair_k) {
                                if (m_packing.is_depleted(pairs[pair_k])) // Note: At least two pairs are depleted.
                                    continue;
                                for (size_t c_i = border[pair_k]; c_i < border[pair_k + 1]; ++c_i) {
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
                                m_packing.subtract_from_potential(candidates[c_i], c_cost);
                            }
                            // Now only a is in the packing.


                            Cost abc_cost = a_cost + b_cost + additional_cost;
                            additional_subgraphs.emplace_back(a_i, a_cost);
                            additional_subgraphs.emplace_back(b_i, b_cost);

                            if (abc_cost > max_cost) {
                                max_cost = abc_cost;
                                max_i.clear();
                                max_i.push_back(std::move(additional_subgraphs));
                            } else if (abc_cost == max_cost) {
                                max_i.push_back(std::move(additional_subgraphs));
                            }
                        }
                    }
                }

                m_packing.add_to_potential(a, a_cost);
                // No candidates remain in the packing.
            }
        }

        assert(!max_i.empty());

        std::uniform_int_distribution<size_t> sample(0, max_i.size() - 1);
        size_t j = sample(m_gen);

        bool changed = max_cost > x_cost || (max_i[j].size() == 1 && candidates[max_i[j].front().first] ==
                                                                     x); // Must be done here, because candidates will be moved into inserted_subgraphs.

        for (auto[y_i, y_cost] : max_i[j]) {
            m_packing.subtract_from_potential(candidates[y_i], y_cost);
            inserted_subgraphs.emplace_back(std::move(candidates[y_i]), y_cost);
        }

        assert(m_packing.is_valid());
        return {changed, max_cost > x_cost};
    }

    /**
     * When an edge is marked, subgraphs in the packing may not have the maximal cost. Also, subgraphs which were
     * previously blocked by a depleted vertex pair, may now be able to be inserted into the packing.
     *
     * Strategy: Remove all subgraphs from the packing incident to the vertex pair. Insert subgraphs to make the packing
     * maximal.
     *
     * @param uv
     */
    void local_fix_marked(VertexPair uv) {

    }

    /**
     * When an edge is edited, subgraphs in the packing may no longer be a forbidden subgraph.
     *
     * Strategy: Remove all subgraphs from the packing incident to the vertex pair. Insert subgraphs to make the packing
     * maximal.
     * @param uv
     */
    void local_fix_edited(VertexPair uv) {

    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKINGLOCALSEARCH_H
