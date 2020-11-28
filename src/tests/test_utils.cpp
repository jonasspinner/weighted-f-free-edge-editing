#include "test_utils.h"


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
        G.set_edge(uv);
    }
    return G;
}

