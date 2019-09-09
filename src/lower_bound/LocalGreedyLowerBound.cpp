//
// Created by jonas on 29.07.19.
//


#include "LocalGreedyLowerBound.h"


namespace LowerBound {
    /**
     * Calculates a lower bound on the costs required to solve the current instance.
     *
     * Inserts forbidden subgraphs into the bound if possible. A subgraph is not inserted if it shares an editable
     * vertex pair with a subgraph already in the bound.
     *
     * @param k Not used
     * @return A lower bound on the costs required to solve the current instance.
     */
    Cost LocalGreedyLowerBound::result(Cost /*k*/) {
        Cost bound_size = 0;
        VertexPairMap<bool> is_in_bound(m_graph.size(), false);

        finder->find([&](Subgraph &&subgraph) {
            bool touches_bound = false;
            for (VertexPair uv : subgraph.vertexPairs()) {
                if (!m_marked[uv] && is_in_bound[uv]) touches_bound = true;
            }
            if (!touches_bound) {
                bound_size += get_subgraph_cost(subgraph, m_marked, m_costs);
                for (VertexPair uv : subgraph.vertexPairs()) {
                    if (!m_marked[uv]) is_in_bound[uv] = true;
                }
            }
            return false;
        });

        return bound_size;
    }
}
