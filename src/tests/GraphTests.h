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
    explicit GraphTests(int seed = 0) : gen(static_cast<unsigned long>(seed)) {}

    template <class It>
    static void check_random_access_iterator(It begin, It end) {
        It i;
        const It j;
        const typename std::iterator_traits<It>::difference_type n;
        static_assert(std::is_same_v<decltype(i += n), It&>);
        static_assert(std::is_same_v<decltype(j +  n), It>);
        static_assert(std::is_same_v<decltype(n +  j), It>);
        static_assert(std::is_same_v<decltype(i -= n), It&>);
        static_assert(std::is_same_v<decltype(j -  n), It>);
        static_assert(std::is_same_v<decltype(j[n]), std::iterator_traits<It>::reference>);

        It a, b;
        (a += n) == b;
        std::addressof(a += n) == std::addressof(b);
        (a + n) == (a += n);
        (a + n) == (n + a);
        a + 0 = a;
        --b == (a + (n - 1));
        ++b;
        (b += -n) = a;
        (b -= n) = a;
        (b - n) == (b -= n);
        a[n] == *b;

        bool(a <= b) == true;
    }

    void vertices_and_for_loop_are_consistent();

    void edges_and_for_loops_are_consistent();

    void vertexPairs_and_for_loops_are_consistent();

    void neighbors_and_for_loops_are_consistent();

    static void iterators_on_empty_Graph_work();

    static void iterators_on_singleton_Graph_work();

    void run();
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_GRAPHTESTS_H
