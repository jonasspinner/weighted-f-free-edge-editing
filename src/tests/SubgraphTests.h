//
// Created by jonas on 09.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHTESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHTESTS_H


#include <random>
#include "../graph/Graph.h"
#include "../graph/Subgraph.h"
#include "test_utils.h"

class SubgraphTests {
    std::mt19937 gen;
public:
    explicit SubgraphTests(int seed = 0) : gen(static_cast<unsigned long>(seed)) {}

    void vertices_and_for_loop_are_consistent() {
        Subgraph subgraph = random_subgraph(10, 20, gen);

        std::vector<Vertex> a, b;

        for (Vertex u : subgraph.vertices())
            a.push_back(u);

        for (size_t i = 0; i < subgraph.size(); ++i)
            b.push_back(subgraph[i]);
    }

    void vertexPairs_and_for_loops_are_consistent() {
        Subgraph subgraph = random_subgraph(10, 20, gen);

        std::vector<VertexPair> a, b;

        for (VertexPair uv : subgraph.vertexPairs())
            a.push_back(uv);

        for (size_t i = 0; i < subgraph.size(); ++i)
            for (size_t j = i + 1; j < subgraph.size(); ++j)
                b.emplace_back(subgraph[i], subgraph[j]);

        expect("vertexPairs() and for loops produce the same output", a, b);
    }

    static void iterators_on_empty_Subgraph_work() {
        Subgraph subgraph{};

        auto vertices = subgraph.vertices();
        expect("empty Subgraph has no vertices", true, vertices.begin() == vertices.end());

        auto vertexPairs = subgraph.vertexPairs();
        expect("empty Subgraph has no vertex pairs", true, vertexPairs.begin() == vertexPairs.end());
    }

    static void iterators_on_singleton_Subgraph_work() {
        Subgraph subgraph{0};

        auto vertices = subgraph.vertices();
        expect("singleton Subgraph has one vertex", true, ++vertices.begin() == vertices.end());

        auto vertexPairs = subgraph.vertexPairs();
        expect("singleton Subgraph has no vertex pairs", true, vertexPairs.begin() == vertexPairs.end());
    }

    void run() {
        std::cout << "\nSubgraphTests"
                     "\n-------------" << std::endl;
        vertices_and_for_loop_are_consistent();
        vertexPairs_and_for_loops_are_consistent();
        iterators_on_empty_Subgraph_work();
        iterators_on_singleton_Subgraph_work();
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHTESTS_H
