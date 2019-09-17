//
// Created by jonas on 03.07.19.
//


#include "Subgraph.h"


/**
 * Returns the minimum cost of all vertex pairs of the subgraph which are not marked. If all vertex pairs are marked
 * invalid_cost is returned. invalid_cost is the maximum of the Cost type.
 *
 * @param subgraph
 * @param marked
 * @param costs
 * @return
 */
Cost get_subgraph_cost(const Subgraph &subgraph, const VertexPairMap<bool> &marked, const VertexPairMap<Cost> &costs) {
    Cost min_cost = invalid_cost;
    for (VertexPair uv : subgraph.vertexPairs())
        if (!marked[uv])
            min_cost = std::min(min_cost, costs[uv]);
    return min_cost;
}
