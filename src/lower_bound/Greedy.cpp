//
// Created by jonas on 29.07.19.
//


#include "Greedy.h"


namespace lower_bound {
    /**
     * Calculates a lower bound on the costs required to solve the current instance.
     *
     * Inserts forbidden subgraphs into the bound if possible. A subgraph is not inserted if it shares an editable
     * vertex pair with a subgraph already in the bound.
     *
     * Complexity:
     *      m = #forbidden subgraphs
     *      O(m) time
     *      O(1) additional space
     *
     * @param k Not used
     * @return A lower bound on the costs required to solve the current instance.
     */
    Cost Greedy::calculate_lower_bound(Cost k) {
        Cost bound_size = 0;
        m_used_in_bound.clear();

        finder->find(m_graph, [&](Subgraph &&subgraph) {

            // Check if the subgraph is adjacent to one already used in the bound.
            bool touches_bound = false;
            for (VertexPair uv : subgraph.vertexPairs())
                if (!m_marked[uv] && m_used_in_bound[uv]) {
                    touches_bound = true;
                    break;
                }

            if (!touches_bound) {
                Cost cost = get_subgraph_cost(subgraph, m_marked, m_costs);
                if (cost == invalid_cost) {
                    bound_size = invalid_cost;
                    return true;
                }
                bound_size += cost;

                for (VertexPair uv : subgraph.vertexPairs())
                    if (!m_marked[uv])
                        m_used_in_bound[uv] = true;
            }
            return bound_size > k;
        });

        return bound_size;
    }
}
