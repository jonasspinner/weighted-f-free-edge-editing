//
// Created by jonas on 29.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_TESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_TESTS_H


#include <iostream>
#include <vector>
#include <algorithm>
#include "../graph/Subgraph.h"

template <typename T>
void expect(const std::string &name, T expected, T actual) {
    if (expected == actual) {
        std::cout << "Test [" << name << "] succeeded\n";
    } else {
        std::cerr << "Test [" << name << "] failed\n\tExpected " << expected << "\n\tGot      " << actual << "\n";
    }
}

template <typename T>
void expect(const std::string &name, std::vector<T> expected, std::vector<T> actual) {
    if (expected == actual) {
        std::cout << "Test [" << name << "] succeeded [" << expected.size() << "]" << std::endl;
    } else {
        std::cerr << "Test [" << name << "] failed\n\tExpected [" << expected.size() << "] {";
        for (const T &t: expected) std::cerr << " " << t;
        std::cerr <<" }\n\tGot      [" << actual.size() << "] {";
        for (const T &t: actual) std::cerr << " " << t;
        std::cerr << " }" << std::endl;
    }
}

template <typename T>
void expect(const std::string &name, std::vector<std::vector<T>> expected, std::vector<std::vector<T>> actual) {
    if (expected == actual) {
        std::cout << "Test [" << name << "] succeeded [" << expected.size() << "]" << std::endl;
    } else {
        std::cerr << "Test [" << name << "] failed\n\tExpected [" << expected.size() << "] {";
        for (auto l : expected)
            for (auto t : l)
                std::cerr << " " << t;
        std::cerr <<" }\n\tGot      [" << actual.size() << "] {";
        for (auto l : actual)
            for (auto t : l)
                std::cerr << " " << t;
        std::cerr << " }" << std::endl;
    }
}


template <typename T>
std::vector<Subgraph> normalize(std::vector<std::vector<T>> list) {
    for (auto& set : list) {
        std::sort(set.begin(), set.end());
    }
    std::sort(list.begin(), list.end());
    return list;
}

std::vector<Subgraph> normalize(std::vector<Subgraph> list) {
    for (auto& subgraph : list) {
        std::sort(subgraph.begin(), subgraph.end());
    }
    std::sort(list.begin(), list.end());
    return list;
}

VertexPair random_vertex_pair(int size, std::mt19937 &gen) {
    std::uniform_int_distribution<Vertex> dist(0, size - 2);
    Vertex u = dist(gen), v = dist(gen);
    v = (u != v) ? v : v + 1;
    return {u, v};
}


Graph random_graph(int size, int n_edges, std::mt19937 &gen) {
    Graph G(size);
    int m = 0;
    while (m < n_edges) {
        auto uv = random_vertex_pair(size, gen);
        G.set_edge(uv);
        m++;
    }
    return G;
}

Subgraph random_subgraph(int size, int graph_size, std::mt19937 &gen) {
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


#endif //WEIGHTED_F_FREE_EDGE_EDITING_TESTS_H
