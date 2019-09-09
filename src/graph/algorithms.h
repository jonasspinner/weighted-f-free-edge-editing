//
// Created by jonas on 08.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ALGORITHMS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ALGORITHMS_H


std::pair<std::vector<Graph>, std::vector<VertexMap<Vertex>>> connectedComponents(const Graph &graph) {
    VertexMap<bool> marked(graph.size(), false);
    std::vector<std::vector<Vertex>> ccVertices;

    for (Vertex x : graph.vertices()) {
        if (marked[x]) continue;

        ccVertices.emplace_back();
        auto &cc = ccVertices.back();

        std::deque<Vertex> stack{x};

        while (!stack.empty()) {
            Vertex u = stack.back();
            stack.pop_back();

            for (Vertex v : graph.neighbors(u))
                if (!marked[v])
                    stack.push_back(v);

            marked[u] = true;
            cc.push_back(u);
        }
    }

    std::vector<Graph> ccs;
    std::vector<VertexMap<Vertex>> mappings;

    for (const auto& vertices : ccVertices) {
        ccs.emplace_back(vertices.size());
        mappings.emplace_back(vertices.size());

        auto &cc = ccs.back();
        auto &m = mappings.back();

        for (Vertex i = 0; i < static_cast<Vertex>(vertices.size()); ++i)
            m[i] = vertices[i];

        for (auto [u, v] : graph.vertexPairs())
            if (graph.hasEdge({m[u], m[v]}))
                cc.setEdge({u, v});
    }

    return {ccs, mappings};
}

bool isConnected(const Graph &graph) {
    VertexMap<bool> marked(graph.size(), false);
    size_t num_marked = 0;

    std::vector<Vertex> stack;
    stack.push_back(0);

    while (!stack.empty()) {
        Vertex u = stack.back();
        stack.pop_back();

        for (Vertex v : graph.neighbors(u))
            if (!marked[v])
                stack.push_back(v);

        marked[u] = true;
        ++num_marked;
    }

    return num_marked == graph.size();
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ALGORITHMS_H
