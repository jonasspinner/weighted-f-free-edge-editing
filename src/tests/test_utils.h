#ifndef WEIGHTED_F_FREE_EDGE_EDITING_TEST_UTILS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_TEST_UTILS_H


#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

#include "../graph/VertexPair.h"
#include "../graph/Graph.h"

template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec) {
    os << "{ ";
    for (const auto &t : vec)
        os << t << " ";
    return os << "}";
}


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
        for (auto l : expected) {
            std::cerr << " {";
            for (auto t : l)
                std::cerr << " " << t;
            std::cerr << " }";
        }
        std::cerr << " }\n\tGot      [" << actual.size() << "] {";
        for (auto l : actual) {
            std::cerr << " {";
            for (auto t : l)
                std::cerr << " " << t;
            std::cerr << " }";
        }
        std::cerr << " }" << std::endl;
    }
}


template <typename T>
std::vector<std::vector<T>> normalize(std::vector<std::vector<T>> list);

VertexPair random_vertex_pair(unsigned size, std::mt19937 &gen);

Graph random_graph(unsigned size, int n_edges, std::mt19937 &gen);

#endif //WEIGHTED_F_FREE_EDGE_EDITING_TEST_UTILS_H
