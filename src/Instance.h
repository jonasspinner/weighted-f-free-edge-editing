//
// Created by jonas on 17.07.19.
//

#ifndef CONCEPT_INSTANCE_H
#define CONCEPT_INSTANCE_H


#include <random>
#include "graph/VertexPairMap.h"

class Instance {
private:
    std::vector<Vertex> permutation;
    std::vector<Vertex> r_permutation;
public:
    Graph graph;
    VertexPairMap<Cost> costs;

    Instance(Graph graph_, VertexPairMap<Cost> costs_) : graph(std::move(graph_)), costs(std::move(costs_)) {}
    Instance(Graph&& graph_, VertexPairMap<Cost>&& costs_) : graph(graph_), costs(costs_) {}

    void apply_permuation(int seed) {
        permutation.resize(graph.size());
        r_permutation.resize(graph.size());

        for (size_t i = 0; i < permutation.size(); ++i)
            permutation[i] = i;

        std::shuffle(permutation.begin(), permutation.end(), std::mt19937_64(seed));

        for (size_t i = 0; i < permutation.size(); ++i)
            r_permutation[permutation[i]] = i;

        const auto &p = permutation;

        Graph new_graph(graph.size());
        VertexPairMap<Cost> new_costs(graph.size());

        for (VertexPair uv : graph.vertexPairs()) {
            VertexPair new_uv{p[uv.u], p[uv.v]};
            if (graph.has_edge(uv))
                new_graph.set_edge(new_uv);
            new_costs[new_uv] = costs[uv];
        }

        graph = std::move(new_graph);
        costs = std::move(new_costs);
    }

    void reverse_permutation() {
        Graph old_graph(graph.size());
        VertexPairMap<Cost> old_costs(graph.size());

        const auto &r_p = r_permutation;

        for (VertexPair uv : graph.vertexPairs()) {
            VertexPair old_uv{r_p[uv.u], r_p[uv.v]};
            if (graph.has_edge(uv))
                old_graph.set_edge(old_uv);
            old_costs[old_uv] = costs[uv];
        }

        graph = std::move(old_graph);
        costs = std::move(old_costs);
    }
};


#endif //CONCEPT_INSTANCE_H
