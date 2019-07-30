//
// Created by jonas on 26.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_BOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_BOUND_H


#include <vector>

#include "graph/Subgraph.h"


class Bound {
    std::vector<Subgraph> subgraphs;
    std::vector<Cost> costs;
    Graph graph;


    [[nodiscard]] const Graph& bound() const {
        return graph;
    }

    void assert_valid() const {
        size_t n_edges = 0;
        for (const auto& subgraph : subgraphs) {
            subgraph.for_all_vertex_pairs([&](VertexPair uv) {
                assert(graph.has_edge(uv));
                n_edges++;
                return false;
            });
        }
        size_t n_edges_graph = 0;
        graph.for_all_vertices([&](Vertex u) {
            n_edges_graph += graph.degree(u);
        });
        n_edges_graph /= 2;
        assert(n_edges == n_edges_graph);
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_BOUND_H
