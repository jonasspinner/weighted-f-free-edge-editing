//
// Created by jonas on 08.08.19.
//


#include "GraphTests.h"

#include <random>

#include "../graph/Graph.h"
#include "test_utils.h"


void GraphTests::vertices_and_for_loop_are_consistent() {
    Graph G = random_graph(10, 40, gen);

    std::vector<Vertex> a, b;

    for (Vertex u : G.vertices())
        a.push_back(u);

    for (Vertex u = 0; u < G.size(); ++u)
        b.push_back(u);

    expect("vertices() and for loop produce the same output", a, b);
}

void GraphTests::edges_and_for_loops_are_consistent() {
    Graph G = random_graph(10, 40, gen);

    std::vector<VertexPair> a, b;

    for (VertexPair uv : G.edges())
        a.push_back(uv);

    for (Vertex u = 0; u < G.size(); ++u)
        for (Vertex v = u + 1; v < G.size(); ++v)
            if (G.has_edge({u, v}))
                b.emplace_back(u, v);

    expect("edges() and for loops produce the same output", a, b);
}

void GraphTests::vertexPairs_and_for_loops_are_consistent() {
    Graph G = random_graph(10, 40, gen);

    std::vector<VertexPair> a, b;

    for (VertexPair uv : G.vertex_pairs())
        a.push_back(uv);

    for (Vertex u = 0; u < G.size(); ++u)
        for (Vertex v = u + 1; v < G.size(); ++v)
            b.emplace_back(u, v);

    expect("vertex_pairs() and for loops produce the same output", a, b);
}

void GraphTests::neighbors_and_for_loops_are_consistent() {
    Graph G = random_graph(10, 40, gen);

    std::vector<std::vector<Vertex>> a(10), b(10);

    for (Vertex u : G.vertices())
        for (Vertex v : G.neighbors(u))
            a[u].push_back(v);

    for (Vertex u = 0; u < G.size(); ++u)
        for (Vertex v = 0; v < G.size(); ++v)
            if (u != v && G.has_edge({u, v}))
                b[u].push_back(v);

    expect("neighbors() and for loops produce the same output", a, b);
}

void GraphTests::iterators_on_empty_Graph_work() {
    Graph graph(0);

    auto vertices = graph.vertices();
    expect("empty Graph has no vertices", true, vertices.begin() == vertices.end());

    auto edges = graph.edges();
    expect("empty Graph has no edges", true, edges.begin() == edges.end());

    auto vertexPairs = graph.vertex_pairs();
    expect("empty Graph has no vertex pairs", true, vertexPairs.begin() == vertexPairs.end());
}

void GraphTests::iterators_on_singleton_Graph_work() {
    Graph graph(1);

    auto vertices = graph.vertices();
    expect("singleton Graph has one vertex", true, ++vertices.begin() == vertices.end());

    auto edges = graph.edges();
    expect("singleton Graph has no edges", true, edges.begin() == edges.end());

    auto vertexPairs = graph.vertex_pairs();
    expect("singleton Graph has no vertex pairs", true, vertexPairs.begin() == vertexPairs.end());
}

void GraphTests::run() {
    std::cout << "\nGraphTests"
                 "\n----------" << std::endl;

    vertices_and_for_loop_are_consistent();
    edges_and_for_loops_are_consistent();
    vertexPairs_and_for_loops_are_consistent();
    neighbors_and_for_loops_are_consistent();
    iterators_on_empty_Graph_work();
    iterators_on_singleton_Graph_work();
}

