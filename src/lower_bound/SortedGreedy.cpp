//
// Created by jonas on 29.07.19.
//


#include "SortedGreedy.h"

#include <vector>


namespace lower_bound {
    /**
     * Calculates a lower bound on the costs required to solve the current instance.
     *
     * Greedily inserts forbidden subgraphs into the bound. Higher minimum editing costs are preferred. A subgraph
     * is not inserted if it shares an editable vertex pair with a subgraph already in the bound.
     *
     * Complexity:
     *      m = #forbidden subgraphs
     *      O(m * log m) time
     *      O(m) additional space
     *
     * @param k Not used
     * @return A lower bound on the costs required to solve the current instance.
     */
    Cost SortedGreedy::calculate_lower_bound(Cost k) {
        Cost max_min_cost = std::numeric_limits<Cost>::min();

        // Find all forbidden subgraphs with editable vertex pairs.
        // The cost for a single forbidden subgraph is the minimum edit cost for an editable vertex pair.
        std::vector<std::pair<Cost, Subgraph>> subgraphs;

        // Only use find and not find_with_duplicates because the first version of the subgraph will be inserted.
        // This relies on the fact that all duplicate version have the same cost. Note: That may change.
        finder->find(m_graph, [&](Subgraph &&subgraph) {
            Cost min_cost = finder->calculate_min_cost(subgraph, m_marked, m_costs);
            subgraphs.emplace_back(min_cost, std::move(subgraph));
            max_min_cost = std::max(max_min_cost, min_cost);
            return max_min_cost > k;
        });

        // If one subgraph has an already large enough cost or if a subgraph is fully marked.
        if (max_min_cost > k)
            return max_min_cost;

        // Sort subgraphs with decreasing costs.
        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

        Cost bound_size = 0;
        m_used_in_bound.clear();

        // Insert forbidden subgraphs with decreasing minimum edit cost into the bound.
        // Only insert a subgraph if it does not share an editable vertex pair with a subgraph already in the bound.
        for (const auto&[cost, subgraph] : subgraphs) {
            if (bound_size > k)
                break;

            // Check if the subgraph is adjacent to one already used in the bound.
            bool touches_bound = finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
                return !m_marked[uv] && m_used_in_bound[uv];
            });

            if (!touches_bound) {
                bound_size += cost;
                finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
                    if (!m_marked[uv])
                        m_used_in_bound[uv] = true;
                    return false;
                });
            }
        }

        return bound_size;
    }

    std::tuple<Cost, VertexPairMap<bool>, std::vector<Subgraph>, std::vector<VertexPair>>
    SortedGreedy::calculate_lower_bound_and_packing() const {

        // Find all forbidden subgraphs with editable vertex pairs.
        // The cost for a single forbidden subgraph is the minimum edit cost for an editable vertex pair.
        std::vector<std::pair<Cost, Subgraph>> subgraphs;
        finder->find(m_graph, [&](Subgraph &&subgraph) {
            Cost min_cost = finder->calculate_min_cost(subgraph, m_marked, m_costs);
            subgraphs.emplace_back(min_cost, std::move(subgraph));
            return false;
        });

        // Sort subgraphs with decreasing costs.
        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

        Cost packing_cost = 0;
        VertexPairMap<bool> covered_by_packing(m_graph.size());
        std::vector<Subgraph> packing;
        std::vector<VertexPair> min_cost_vertex_pairs;

        // Insert forbidden subgraphs with decreasing minimum edit cost into the bound.
        // Only insert a subgraph if it does not share an editable vertex pair with a subgraph already in the bound.
        for (const auto&[cost, subgraph] : subgraphs) {

            // Check if the subgraph is adjacent to one already used in the bound.
            bool touches_bound = false;
            for (VertexPair uv : subgraph.vertexPairs())
                if (!m_marked[uv] && covered_by_packing[uv]) {
                    touches_bound = true;
                    break;
                }

            if (!touches_bound) {
                VertexPair xy{-0, static_cast<Vertex>(-1)};
                Cost min_cost = invalid_cost;
                for (VertexPair uv : subgraph.vertexPairs()) {
                    if (!m_marked[uv]) {
                        covered_by_packing[uv] = true;
                        if (m_costs[uv] < min_cost) {
                            min_cost = m_costs[uv];
                            xy = uv;
                        }
                    }
                }

                packing_cost += cost;
                packing.push_back(subgraph);
                min_cost_vertex_pairs.push_back(xy);
            }
        }

        return {packing_cost, covered_by_packing, packing, min_cost_vertex_pairs};
    }
}
