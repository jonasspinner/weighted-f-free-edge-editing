//
// Created by jonas on 08.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GRAPHTESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GRAPHTESTS_H


#include <random>
#include "../graph/Graph.h"
#include "Tests.h"

class GraphTests {
    std::mt19937 gen;
public:
    /**
     * Tests for checking consistency between range iteration and templated loops with lambda functions and tests for edge cases for graphs.
     * @param seed
     */
    explicit GraphTests(int seed = 0) : gen(seed) {}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    void vertices_and_for_all_vertices_are_consistent() {
        Graph G = random_graph(10, 40, gen);

        std::vector<Vertex> a, b;

        for (Vertex u : G.vertices()) a.push_back(u);

        G.for_all_vertices([&](Vertex u) {
            b.push_back(u);
            return false;
        });

        expect("vertices() and for_all_vertices() produce the same output", a, b);
    }

    void edges_and_for_all_edges_are_consistent() {
        Graph G = random_graph(10, 40, gen);

        std::vector<VertexPair> a, b;

        for (VertexPair uv : G.edges()) a.push_back(uv);

        G.for_all_edges([&](VertexPair uv) {
            b.push_back(uv);
            return false;
        });

        expect("edges() and for_all_edges() produce the same output", a, b);
    }

    void vertexPairs_and_for_all_vertex_pairs_are_consistent() {
        Graph G = random_graph(10, 40, gen);

        std::vector<VertexPair> a, b;

        for (VertexPair uv : G.vertexPairs()) a.push_back(uv);

        G.for_all_vertex_pairs([&](VertexPair uv) {
            b.push_back(uv);
            return false;
        });

        expect("vertexPairs() and for_all_vertex_pairs() produce the same output", a, b);
    }

    void neighbors_and_for_all_neighbors_are_consistent() {
        Graph G = random_graph(10, 40, gen);

        std::vector<std::vector<Vertex>> a(10), b(10);

        for (Vertex u : G.vertices())
            for (Vertex v : G.neighbors(u))
                a[u].push_back(v);

        G.for_all_vertices([&](Vertex u) {
            G.for_neighbors_of(u, [&](Vertex v) {
                b[u].push_back(v);
                return false;
            });
            return false;
        });

        expect("neighbors() and for_all_neighbors() produce the same output", a, b);
    }
#pragma GCC diagnostic pop

    static void iterators_on_empty_Graph_work() {
        Graph graph(0);

        auto vertices = graph.vertices();
        expect("empty Graph has no vertices", true, vertices.begin() == vertices.end());

        auto edges = graph.edges();
        expect("empty Graph has no edges", true, edges.begin() == edges.end());

        auto vertexPairs = graph.vertexPairs();
        expect("empty Graph has no vertex pairs", true, vertexPairs.begin() == vertexPairs.end());
    }

    static void iterators_on_singleton_Graph_work() {
        Graph graph(1);

        auto vertices = graph.vertices();
        expect("singleton Graph has one vertex", true, ++vertices.begin() == vertices.end());

        auto edges = graph.edges();
        expect("singleton Graph has no edges", true, edges.begin() == edges.end());

        auto vertexPairs = graph.vertexPairs();
        expect("singleton Graph has no vertex pairs", true, vertexPairs.begin() == vertexPairs.end());
    }

    void run() {
        std::cout << "\nGraphTests"
                     "\n----------" << std::endl;

        vertices_and_for_all_vertices_are_consistent();
        edges_and_for_all_edges_are_consistent();
        vertexPairs_and_for_all_vertex_pairs_are_consistent();
        neighbors_and_for_all_neighbors_are_consistent();
        iterators_on_empty_Graph_work();
        iterators_on_singleton_Graph_work();
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_GRAPHTESTS_H
