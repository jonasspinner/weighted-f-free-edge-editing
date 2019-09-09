//
// Created by jonas on 22.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SYNTHETIC_GRAPHS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SYNTHETIC_GRAPHS_H


#include "VertexMap.h"


/**
 * Create
 * @return
 */
Instance create_synthetic_cluster_graph(size_t num_nodes, int seed) {

    std::vector<std::vector<Vertex>> clusters;

    std::mt19937 gen(seed);

    size_t n = num_nodes;
    while (n > 0) {
        std::uniform_int_distribution<Vertex> dist(1, n);
        size_t cluster_size = dist(gen);

        std::cout << "cluster size " << cluster_size << "\n";
        n -= cluster_size;
        std::cout << "n' " << n << "\n";
        clusters.emplace_back();
        for (size_t i = 0; i < cluster_size; ++i)
            clusters.back().push_back(n + i);
    }


    std::normal_distribution<double> intra_cluster_dist(100, 80);
    std::normal_distribution<double> inter_cluster_dist(-100, 80);

    VertexPairMap<double> weights(num_nodes);

    for (size_t c = 0; c < clusters.size(); ++c) {
        // intra cluster costs
        for (size_t j = 0; j < clusters[c].size(); ++j)
            for (size_t k = j + 1; k < clusters[c].size(); ++k) {
                Vertex u = clusters[c][j];
                Vertex v = clusters[c][k];
                weights[{u, v}] = intra_cluster_dist(gen);
            }

        // inter cluster costs
        for (size_t c_2 = c + 1; c_2 < clusters.size(); ++c_2)
            for (Vertex u : clusters[c])
                for (Vertex v : clusters[c_2])
                    weights[{u, v}] = inter_cluster_dist(gen);
    }


    Graph graph(num_nodes);
    VertexPairMap<Cost> costs(num_nodes);
    for (VertexPair uv : graph.vertexPairs()) {
        if (weights[uv] >= 0) graph.setEdge(uv);
        costs[uv] = std::abs(weights[uv]);
    }

    std::cout << graph;
    std::cout << costs;

    std::stringstream name;
    name << "synthetic-cluster-graph-size-" << num_nodes;
    return {name.str(), graph, costs};
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SYNTHETIC_GRAPHS_H
