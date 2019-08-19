//
// Created by jonas on 17.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GRAPHSTATISTICS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GRAPHSTATISTICS_H


#include <cstddef>
#include "graph/Graph.h"
#include "graph/Subgraph.h"
#include "finder/NaiveC4P4.h"

class GraphStatistics {
    size_t n_vertices;
    size_t n_edges;


    static GraphStatistics calculate(const Graph &graph) {
        std::vector<size_t> degrees(graph.size());
        for (Vertex u : graph.vertices())
            degrees[u] = graph.degree(u);

        std::vector<Subgraph> subgraphs;
        Finder::NaiveC4P4 finder(graph);
        finder.find([&](Subgraph &&subgraph) {
            subgraphs.push_back(std::move(subgraph));
            return false;
        });

        Graph SG_graph(subgraphs.size());
        for (Vertex i = 0; i < subgraphs.size(); ++i) {
            auto sg_i_vp = subgraphs[i].vertexPairs();
            for (Vertex j = i + 1; j < subgraphs.size(); ++j) {
                bool adjacent = std::any_of(sg_i_vp.begin(), sg_i_vp.end(), [&](VertexPair uv) {
                    for (VertexPair xy : subgraphs[j].vertexPairs())
                        if (xy == uv) return true;
                    return false;
                });
                if (adjacent)
                    SG_graph.set_edge({i, j});
            }
        }

        std::vector<size_t> sg_degrees(SG_graph.size());
        for (Vertex i : SG_graph.vertices())
            sg_degrees[i] = SG_graph.degree(i);
    }

    size_t complexity() {
        return n_vertices * n_edges;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_GRAPHSTATISTICS_H
