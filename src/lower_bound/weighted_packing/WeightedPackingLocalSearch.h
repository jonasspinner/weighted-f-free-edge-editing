#ifndef WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKINGLOCALSEARCH_H
#define WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKINGLOCALSEARCH_H

#include "WeightedPacking.h"

#include "robin_hood.h"
#include <unordered_set>


class WeightedPackingLocalSearch : public LowerBoundI {
    class State {
    public:
        robin_hood::unordered_map<Subgraph, Cost> m_subgraphs_in_packing;

        [[nodiscard]] size_t num_subgraphs_in_packing() const {
            return m_subgraphs_in_packing.size();
        }

        [[nodiscard]] const auto &subgraphs() const {
            return m_subgraphs_in_packing;
        }
    };

    const Graph &m_graph;
    const VertexPairMap<Cost> &m_costs;
    const VertexPairMap<bool> &m_marked;

    WeightedPacking m_packing;
    std::vector<std::unique_ptr<State>> m_states;

public:
    WeightedPackingLocalSearch(const Instance &instance, const VertexPairMap<bool> &marked,
                               const SubgraphStats &subgraph_stats, std::shared_ptr<FinderI> finder_ref) :
            LowerBoundI(std::move(finder_ref)), m_graph(instance.graph), m_costs(instance.costs), m_marked(marked),
            m_packing(instance, marked, subgraph_stats) {
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
        greedy_initialize(state, k);
        local_search(state);
        return m_packing.cost();
    }

    void local_search(State &state) {
        size_t max_iter = 5;
        size_t max_rounds_no_improvement = 5;
        size_t num_rounds_no_improvement = 0;
        std::vector<std::pair<Subgraph, Cost>> total_removed, total_inserted;

        for (size_t iter = 0; iter < max_iter && num_rounds_no_improvement < max_rounds_no_improvement; ++iter) {
            bool found_improvement = false;
            total_removed.clear();
            total_inserted.clear();

            for (const auto &[subgraph, cost] : state.subgraphs()) {
                auto[changed, improved, removed_from_packing, inserted_into_packing] = \
                    find_one_two_improvement(state, subgraph, cost);
                found_improvement |= improved;

                std::move(removed_from_packing.begin(), removed_from_packing.end(),
                          std::back_inserter(total_removed));
                std::move(inserted_into_packing.begin(), inserted_into_packing.end(),
                          std::back_inserter(total_inserted));
            }

            // Update state after iteration. During a iteration, the state is stale and updates are only represented in
            // m_packing.
            // Insert before remove, because both could cancel each other.
            for (auto&&[subgraph, cost] : total_inserted) {
                // If the subgraph is not yet in the map, the initial value is 0.
                state.m_subgraphs_in_packing[std::move(subgraph)] += cost;
            }
            for (auto&&[subgraph, cost] : total_removed) {
                // If the subgraph is fully removed, i.e. its cost is zero, it is removed from the map.
                auto it = state.m_subgraphs_in_packing.find(subgraph);
                assert(it->second >= cost);
                it->second -= cost;
                if (it->second == 0) {
                    state.m_subgraphs_in_packing.erase(it);
                }
            }

            num_rounds_no_improvement = found_improvement ? 0 : num_rounds_no_improvement + 1;
        }
    }

    bool greedy_initialize(State &state, Cost k) {
        std::vector<std::pair<Subgraph, Cost>> subgraph_heap; // Note: Could turn this into an attribute and call clear() at beginning. This would prevent reallocating memory.
        Cost lower_bound = 0;
        Cost max_min_cost = std::numeric_limits<Cost>::min();

        bool unsolvable = finder->find_with_duplicates(m_graph, [&](Subgraph &&subgraph) {
            Cost initial_min_cost = m_packing.calculate_min_cost(subgraph);

            subgraph_heap.emplace_back(std::move(subgraph), initial_min_cost);

            max_min_cost = std::max(max_min_cost, initial_min_cost);
            return max_min_cost > k;
        });

        if (unsolvable)
            return false;

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
                state.m_subgraphs_in_packing.emplace(std::move(subgraph_heap.back()));

                subgraph_heap.pop_back();
            } else {
                // The editing cost is no longer up to date. The subgraph may be inserted in the future.
                subgraph_heap.back().second = current_min_cost;
                std::push_heap(subgraph_heap.begin(), subgraph_heap.end());
            }
        }

        assert(m_packing.is_maximal());
        return true;
    }


    std::tuple<bool, bool, std::vector<std::pair<Subgraph, Cost>>, std::vector<std::pair<Subgraph, Cost>>>
    find_one_two_improvement(State &state, const Subgraph &x, Cost x_cost) {
        std::vector<std::pair<Subgraph, Cost>> removed_from_packing, inserted_into_packing; // Note: could take total_removed and total_inserted as mutable reference parameters.

        m_packing.add_to_potential(x, x_cost);
        removed_from_packing.emplace_back(x, x_cost);

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

                if (a_cost > max_cost) {
                    max_cost = a_cost;
                    max_i = {{{a_i, a_cost}}};
                } else if (a_cost == max_cost) {
                    max_i.push_back({{a_i, a_cost}});
                }

                for (size_t pair_j = pair_i + 1; pair_j < pairs.size(); ++pair_j) {
                    if (m_packing.is_depleted(pairs[pair_j]))
                        continue;
                    for (size_t b_i = border[pair_j]; b_i < border[pair_j + 1]; ++b_i) {
                        const Subgraph &b = candidates[b_i];

                        auto b_cost = m_packing.calculate_min_cost(b);

                        if (b_cost > 0) {
                            // b can be inserted.

                            if (a_cost + b_cost > max_cost) {
                                max_cost = a_cost + b_cost;
                                max_i = {{{a_i, a_cost}, {b_i, b_cost}}};
                            } else if (a_cost + b_cost == max_cost) {
                                max_i.push_back({{a_i, a_cost},
                                                 {b_i, b_cost}});
                            }
                        }
                    }
                }

                m_packing.add_to_potential(a, a_cost);
            }
        }

        if (max_i.empty()) {
            // No improvement has been found. x is still optimal.
            m_packing.subtract_from_potential(x, x_cost);
            // TODO: Should not happen, as x is also a candidate.
            assert(false);
        } else {
            // Plateau search
            size_t j = 0; // currently first element.
            // TODO: Plateau search strategy.
            for (auto[y_i, y_cost] : max_i[j]) {
                m_packing.subtract_from_potential(candidates[y_i], y_cost);
                inserted_into_packing.emplace_back(std::move(candidates[y_i]), y_cost);
            }
        }

        assert(m_packing.is_valid());
        return {max_cost > x_cost || removed_from_packing != inserted_into_packing, max_cost > x_cost,
                removed_from_packing, inserted_into_packing};
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
