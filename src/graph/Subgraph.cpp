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

void verify(const Subgraph &subgraph, const Graph &graph) {
    assert(subgraph.m_vertices.size() == 4);

    const Subgraph &S = subgraph;

    assert(S[0] != S[1]);
    assert(S[0] != S[2]);
    assert(S[0] != S[3]);
    assert(S[1] != S[2]);
    assert(S[1] != S[3]);
    assert(S[2] != S[3]);

    bool is_one = false;
    std::vector<size_t> map = {0, 1, 2, 3};
    do {
        size_t i = map[0], j = map[1], k = map[2], l = map[3];
        bool is_path =
                graph.hasEdge({S[i], S[j]}) &&
                !graph.hasEdge({S[i], S[k]}) &&
                !graph.hasEdge({S[i], S[l]}) &&
                graph.hasEdge({S[j], S[k]}) &&
                !graph.hasEdge({S[j], S[l]}) &&
                graph.hasEdge({S[k], S[l]});
        bool is_cycle =
                graph.hasEdge({S[i], S[j]}) &&
                !graph.hasEdge({S[i], S[k]}) &&
                graph.hasEdge({S[i], S[l]}) &&
                graph.hasEdge({S[j], S[k]}) &&
                !graph.hasEdge({S[j], S[l]}) &&
                graph.hasEdge({S[k], S[l]});
        is_one |= is_path || is_cycle;
    } while (std::next_permutation(map.begin(), map.end()));

    if (!is_one) {
        std::cout << subgraph << " E = {";
        for (VertexPair uv : subgraph.vertexPairs()) {
            if (graph.hasEdge(uv)) std::cout << " " << uv;
        }
        std::cout << " }\n";
    }
    assert(is_one);
}