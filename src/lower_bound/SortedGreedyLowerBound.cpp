//
// Created by jonas on 29.07.19.
//


#include "SortedGreedyLowerBound.h"

#include <vector>


namespace LowerBound {
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
    Cost SortedGreedyLowerBound::calculate_lower_bound(Cost k) {
        Cost max_min_cost = std::numeric_limits<Cost>::min();

        // Find all forbidden subgraphs with editable vertex pairs.
        // The cost for a single forbidden subgraph is the minimum edit cost for an editable vertex pair.
        std::vector<std::pair<Cost, Subgraph>> subgraphs;
        finder->find([&](Subgraph &&subgraph) {
            Cost min_cost = get_subgraph_cost(subgraph, m_marked, m_costs);
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
            bool touches_bound = false;
            for (VertexPair uv : subgraph.vertexPairs())
                if (!m_marked[uv] && m_used_in_bound[uv]) {
                    touches_bound = true;
                    break;
                }

            if (!touches_bound) {
                bound_size += cost;
                for (VertexPair uv : subgraph.vertexPairs())
                    if (!m_marked[uv])
                        m_used_in_bound[uv] = true;
            }
        }

        return bound_size;
    }
}
