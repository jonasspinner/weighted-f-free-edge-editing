//
// Created by jonas on 29.07.19.
//


#include "GlobalGreedyLowerBound.h"

#include <vector>


namespace LowerBound {
    /**
     * Calculates a lower bound on the costs required to solve the current instance.
     *
     * Greedily inserts forbidden subgraphs into the bound. Higher minimum editing costs are preferred. A subgraph
     * is not inserted if it shares an editable vertex pair with a subgraph already in the bound.
     *
     * @param k Not used
     * @return A lower bound on the costs required to solve the current instance.
     */
    Cost GlobalGreedyLowerBound::result(Cost /*k*/) {
        // Find all forbidden subgraphs with editable vertex pairs
        // The cost for a single forbidden subgraph is the minimum edit cost for an editable vertex pair
        std::vector<std::pair<Cost, Subgraph>> subgraphs;
        finder->find([&](Subgraph &&subgraph) {
            Cost min_cost = get_subgraph_cost(subgraph, m_marked, m_costs);
            subgraphs.emplace_back(min_cost, std::move(subgraph));
            return false;
        });

        // Sort subgraphs with decreasing costs
        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

        if (!subgraphs.empty() && subgraphs[0].first == invalid_cost) return 0;

        Cost bound_size = 0;
        VertexPairMap<bool> is_in_bound(m_graph.size(), false);

        // Insert forbidden subgraphs with decreasing minimum edit cost into the bound
        // Only insert a subgraph if it does not share an editable vertex pair with a subgraph already in the bound
        for (const auto&[cost, subgraph] : subgraphs) {
            bool touches_bound = false;
            for (VertexPair uv : subgraph.vertexPairs()) {
                if (!m_marked[uv] && is_in_bound[uv])
                    touches_bound = true;
            }

            if (!touches_bound) {
                bound_size += cost;
                for (VertexPair uv : subgraph.vertexPairs()) {
                    if (!m_marked[uv])
                        is_in_bound[uv] = true;
                }
            }
        }

        return bound_size;
    }
}
