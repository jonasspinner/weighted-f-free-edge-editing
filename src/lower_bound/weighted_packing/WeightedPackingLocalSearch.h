#ifndef WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKINGLOCALSEARCH_H
#define WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKINGLOCALSEARCH_H

#include "WeightedPacking.h"


class WeightedPackingLocalSearch : public LowerBoundI {

    const Graph &m_graph;
    const VertexPairMap<bool> &m_marked;
    const VertexPairMap<Cost> &m_costs;

    class State {
    public:
        WeightedPacking m_packing;
        std::vector<std::pair<Subgraph, Cost>> m_subgraphs_in_packing;

        [[nodiscard]] size_t num_subgraphs_in_packing() const {
            return m_subgraphs_in_packing.size();
        }

        [[nodiscard]] const auto &subgraph_in_packing(size_t index) const {
            return m_subgraphs_in_packing[index];
        }

        void remove(size_t index) {
            std::swap(m_subgraphs_in_packing[index], m_subgraphs_in_packing.back());
            const auto&[subgraph, cost] = m_subgraphs_in_packing.back();
            m_packing.add_to_potential(subgraph, cost);
            m_subgraphs_in_packing.pop_back();
        }

        void replace(size_t index, Subgraph &&subgraph, Cost cost) {
            const auto&[old_subgraph, old_cost] = m_subgraphs_in_packing[index];
            m_packing.add_to_potential(old_subgraph, old_cost);
            m_packing.subtract_from_potential(subgraph, cost);
            m_subgraphs_in_packing[index] = {std::move(subgraph), cost};
        }

        void add(Subgraph &&subgraph, Cost cost) {
            m_packing.subtract_from_potential(subgraph, cost);
            m_subgraphs_in_packing.emplace_back(std::move(subgraph), cost);
        }

        void add(std::pair<Subgraph, Cost> &&element) {
            const auto&[subgraph, cost] = element;
            m_packing.subtract_from_potential(subgraph, cost);
            m_subgraphs_in_packing.push_back(std::move(element));
        }

        [[nodiscard]] bool is_valid() const {
            return m_packing.is_valid();
        }
    };

    std::vector<std::unique_ptr<State>> m_states;

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
        return state.m_packing.cost();
    }

    void local_search(State &state) {
        size_t max_iter = 5;
        size_t max_rounds_no_improvement = 5;
        size_t num_rounds_no_improvement = 0;

        for (size_t iter = 0; iter < max_iter && num_rounds_no_improvement < max_rounds_no_improvement; ++iter) {
            bool found_improvement = false;
            for (size_t index = 0; index < state.m_subgraphs_in_packing.size(); ++index) {
                found_improvement |= find_one_two_improvement(state, index);
            }
            num_rounds_no_improvement = found_improvement ? 0 : num_rounds_no_improvement + 1;
        }
    }

    bool greedy_initialize(State &state, Cost k) {
        auto &packing = state.m_packing;

        std::vector<std::pair<Subgraph, Cost>> subgraph_heap;
        Cost lower_bound = 0;
        Cost max_min_cost = std::numeric_limits<Cost>::min();

        bool unsolvable = finder->find_with_duplicates(m_graph, [&](Subgraph &&subgraph) {
            Cost initial_min_cost = packing.calculate_min_cost(subgraph);

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

            Cost current_min_cost = packing.calculate_min_cost(subgraph);

            if (current_min_cost == 0) {
                // The subgraph cannot be inserted.
                subgraph_heap.pop_back();
            } else if (current_min_cost == cost) {
                // The editing cost is maximal from all remaining subgraphs in the queue. Increase the lower bound
                // and update the remaining cost matrix.
                lower_bound += current_min_cost;

                state.add(std::move(subgraph_heap.back()));
                subgraph_heap.pop_back();
            } else {
                // The editing cost is no longer up to date. The subgraph may be inserted in the future.
                subgraph_heap.back().second = current_min_cost;
                std::push_heap(subgraph_heap.begin(), subgraph_heap.end());
            }
        }

        assert(is_maximal());
        return true;
    }


    bool find_one_two_improvement(State &state, size_t index) {
        auto &packing = state.m_packing;

        const auto &[x, x_cost] = state.subgraph_in_packing(index);

        packing.add_to_potential(x, x_cost);


        const auto pairs = packing.get_neighbor_pairs(x);
        auto[candidates, border] = packing.get_neighbors(pairs);

        // The information about the best solution.
        Cost max_cost = x_cost;
        std::vector<std::pair<size_t, Cost>> max_i;


        for (size_t pair_i = 0; pair_i < pairs.size(); ++pair_i) {
            for (size_t a_i = border[pair_i]; a_i < border[pair_i + 1]; ++a_i) {
                const Subgraph &a = candidates[a_i];

                auto a_cost = packing.calculate_min_cost(a);
                packing.subtract_from_potential(a, a_cost);

                if (a_cost > max_cost) {
                    max_cost = a_cost;
                    max_i = {{a_i, a_cost}};
                }

                for (size_t pair_j = pair_i + 1; pair_j < pairs.size(); ++pair_j) {
                    if (packing.is_depleted(pairs[pair_j]))
                        continue;
                    for (size_t b_i = border[pair_j]; b_i < border[pair_j + 1]; ++b_i) {
                        const Subgraph &b = candidates[b_i];

                        auto b_cost = packing.calculate_min_cost(b);

                        if (b_cost > 0) {
                            // b can be inserted.

                            if (a_cost + b_cost > max_cost) {
                                max_cost = a_cost + b_cost;
                                max_i = {{a_i, a_cost},
                                         {b_i, b_cost}};
                            }
                        }
                    }
                }

                packing.add_to_potential(a, a_cost);
            }
        }

        if (max_i.empty()) {
            // No improvement has been found. x is still optimal.
            packing.subtract_from_potential(x, x_cost);
        } else if (max_cost == x_cost) {
            // Plateau search
        } else {
            for (auto[y_i, y_cost] : max_i)
                packing.subtract_from_potential(candidates[y_i], y_cost);
        }

        assert(packing.is_valid());
        return max_cost > x_cost;
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
