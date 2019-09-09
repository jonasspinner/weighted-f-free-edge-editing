//
// Created by jonas on 29.07.19.
//


#include "test_utils.h"


std::vector<Subgraph> normalize(std::vector<Subgraph> list) {
    for (auto& subgraph : list) {
        subgraph.sortVertices();
    }
    std::sort(list.begin(), list.end());
    return list;
}

VertexPair random_vertex_pair(unsigned size, std::mt19937 &gen) {
    std::uniform_int_distribution<Vertex> dist(0, size - 2);
    Vertex u = dist(gen), v = dist(gen);
    v = (u != v) ? v : v + 1;
    return {u, v};
}


Graph random_graph(unsigned size, int n_edges, std::mt19937 &gen) {
    Graph G(size);
    for (int m = 0; m < n_edges; ++m) {
        auto uv = random_vertex_pair(size, gen);
        G.setEdge(uv);
    }
    return G;
}

Subgraph random_subgraph(unsigned size, unsigned graph_size, std::mt19937 &gen) {
    std::uniform_int_distribution<Vertex> dist(0, graph_size - 1);
    std::vector<bool> marked(graph_size);

    Subgraph subgraph{};
    while (subgraph.size() < size) {
        Vertex u = dist(gen);
        while (marked[u]) u = u + 1 % graph_size;
        subgraph.push_back(u);
        marked[u] = true;
    }
    return subgraph;
}
