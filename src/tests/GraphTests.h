//
// Created by jonas on 08.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GRAPHTESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GRAPHTESTS_H


#include <random>
#include "../graph/Graph.h"
#include "test_utils.h"

class GraphTests {
    std::mt19937 gen;
public:
    /**
     * Tests for checking consistency between range iteration and templated loops with lambda functions and tests for edge cases for graphs.
     * @param seed
     */
    explicit GraphTests(int seed = 0) : gen(seed) {}

    void vertices_and_for_all_vertices_are_consistent();

    void edges_and_for_all_edges_are_consistent();

    void vertexPairs_and_for_all_vertex_pairs_are_consistent();

    void neighbors_and_for_all_neighbors_are_consistent();

    static void iterators_on_empty_Graph_work();

    static void iterators_on_singleton_Graph_work();

    void run();
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_GRAPHTESTS_H
