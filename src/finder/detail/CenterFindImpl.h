//
// Created by jonas on 09.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_CENTERFINDIMPL_H
#define WEIGHTED_F_FREE_EDGE_EDITING_CENTERFINDIMPL_H


namespace detail {
    /**
     * Implementation for finding C_k, P_k induced subgraphs.
     *
     * Complexity:
     *      O(n^{k-1} + p_k(G) + k * c_k(G)) to list all P_k's of a graph G, where k >= 4
     *      O(n + m + p_3(G) + co-p_3(G)) to list all P_3's of a graph G
     *
     * Reference:
     *      Chính T. Hoàng, Marcin Kamiński, Joe Sawada, R. Sritharan, Finding and listing induced paths and cycles,
     *      https://doi.org/10.1016/j.dam.2012.01.024.
     *
     * @tparam k size of the subgraphs C_k, P_k
     * @tparam with_cycles whether to search for C_k
     */
    template <int k, bool with_cycles>
    class CenterFindImpl {
    public:
        static_assert(k >= 4);

        /**
         * Find paths (and cycles) with length of at least 4. Each path (and cycle) is listed exactly once.
         *
         * The algorithm recurses on paths of size k - 4 for k >= 6 and on k - 2 for k >= 4.
         *
         * The callback is called with (Subgraph&&, Vertex).
         *
         * @tparam k Length of the paths (or cycles).
         * @tparam with_cycles Whether to search for cycles.
         * @param graph
         * @return
         */
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            Graph::AdjRow A(graph.size()), B(graph.size());

            if constexpr (k >= 6) {
                Graph::AdjRow C(graph.size()), A_p(graph.size()), B_p(graph.size());

                /** Recursion on finding P_{k-4} **/
                return CenterFindImpl<k - 4, false>::find(graph, [&](Subgraph&& P, Vertex P_min_vertex) {
                    // P = p_1, ..., p_{k-4}
                    assert(P.size() == k - 4);


                    C = graph.all_vertices(); // non adjacent to P - {p_1, p_{k-4}} = p_2, ..., p_{k-5}
                    for (unsigned i = 1; i < k - 5; ++i) { C &= non_neighbors(P[i]); }

                    A = neighbors(P[0]); // adjacent to p_1 and non adjacent to P - {p_1} = p_2, ..., p_{k-4}
                    A &= C; A &= non_neighbors(P[k-5]);

                    B = neighbors(P[k - 5]); // adjacent to p_{k-4} and non adjacent to P - {p_{k-4}} = p_1, ..., p_{k-5}
                    B &= C; B &= non_neighbors(P[0]);

                    // non adjacent to P = p_1, ..., p_{k-4}
                    C &= non_neighbors(P[0]);
                    C &= non_neighbors(P[k-5]);


#ifndef NDEBUG
                    for (Vertex a : Graph::iterate(A)) {
                        assert(valid_edge({a, P[0]}));
                        for (unsigned i = 1; i < k-4; ++i)
                            assert(valid_non_edge({a, P[i]}));
                    }
                    for (Vertex b : Graph::iterate(B)) {
                        assert(valid_edge({b, P[k-5]}));
                        for (unsigned i = 0; i < k-5; ++i)
                            assert(valid_non_edge({b, P[i]}));
                    }
                    for (Vertex c : Graph::iterate(C)) {
                        for (Vertex p : P.vertices())
                            assert(valid_non_edge({c, p}));
                    }
#endif

                    // for each (a, b) \in AxB
                    for (Vertex a : Graph::iterate(A)) {
                        for (Vertex b : Graph::iterate(B)) {
#ifndef NDEBUG
                            assert(a != b);
                            for (Vertex p : P.vertices()) { assert(a != p); assert(b != p); }
#endif
                            if (valid_non_edge({a, b})) {

                                A_p = C; A_p &= neighbors(a); A_p &= non_neighbors(b); // subset of C adjacent to a but not b
                                B_p = C; B_p &= neighbors(b); B_p &= non_neighbors(a); // subset of C adjacent to b but not a

                                // for each (u, v) \in A'xB'
                                for (Vertex u : Graph::iterate(A_p)) {
                                    for (Vertex v : Graph::iterate(B_p)) {
                                        Vertex min_vertex = std::min({P_min_vertex, a, b, u, v});
#ifndef NDEBUG
                                        assert(u != v); assert(u != a); assert(u != b); assert(v != a); assert(v != b);
                                        for (Vertex p : P.vertices()) { assert(u != p); assert(v != p); }
#endif
                                        if (valid_non_edge({u, v})) {
                                            // P' = uaPbv
                                            Subgraph P_p({u, a}, P, {b, v});

#ifndef NDEBUG
                                            // assert that P' is a path of length k
                                            assert(P_p.size() == k);
                                            for (unsigned i = 0; i < k; ++i)
                                                for (unsigned j = i + 1; j < k; ++j)
                                                    if (j - i == 1) {
                                                        assert(valid_edge({P_p[i], P_p[j]}));
                                                    } else {
                                                        assert(valid_non_edge({P_p[i], P_p[j]}));
                                                    }
#endif

                                            if (callback(std::move(P_p), min_vertex)) return true;

                                        } else if (valid_edge({u, v}) && with_cycles) { // See 3.4 Listing C_k
                                            // P' = Pbvua
                                            Subgraph P_p({}, P, {b, v, u, a});

                                            if (P_p[0] == min_vertex && P_p[1] < P_p[k-1]) {

#ifndef NDEBUG
                                                // assert that P' is a cycle of length k
                                                assert(P_p.size() == k);
                                                for (unsigned i = 0; i < k; ++i)
                                                    for (unsigned j = i + 1; j < k; ++j)
                                                        if (j - i == 1 || j - i == k-1) {
                                                            assert(valid_edge({P_p[i], P_p[j]}));
                                                        } else {
                                                            assert(valid_non_edge({P_p[i], P_p[j]}));
                                                        }
#endif

                                                if (callback(std::move(P_p), min_vertex)) return true;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // no early exit
                    return false;
                }, neighbors, non_neighbors, valid_edge, valid_non_edge);
            } else if constexpr (k >= 4) {

                /** Recursion on finding P_{k-2} **/
                return CenterFindImpl<k - 2, false>::find(graph, [&](Subgraph &&P, Vertex P_min_vertex) {
                    // P = p_1, ..., p_{k-2}
                    assert(P.size() == k - 2);

                    A = neighbors(P[0]); // adjacent to p_1 and non adjacent to P - {p_1}
                    for (unsigned i = 1; i < k - 2; ++i) { A &= non_neighbors(P[i]); }

                    B = neighbors(P[k - 3]); // adjacent to p_{k-2} and non adjacent to P - {p_{k-2}}
                    for (unsigned i = 0; i < k - 3; ++i) { B &= non_neighbors(P[i]); }


#ifndef NDEBUG
                    for (Vertex a : Graph::iterate(A)) {
                        assert(valid_edge({a, P[0]}));
                        for (unsigned i = 1; i < k-2; ++i)
                            assert(valid_non_edge({a, P[i]}));
                    }
                    for (Vertex b : Graph::iterate(B)) {
                        assert(valid_edge({b, P[k-3]}));
                        for (unsigned i = 0; i < k-3; ++i)
                            assert(valid_non_edge({b, P[i]}));
                    }
#endif

                    // for each (a, b) \in AxB
                    for (Vertex a : Graph::iterate(A)) {
                        for (Vertex b : Graph::iterate(B)) {
                            Vertex min_vertex = std::min({P_min_vertex, a, b});
#ifndef NDEBUG
                            assert(a != b);
                            for (Vertex p : P.vertices()) { assert(a != p); assert(b != p); }
#endif
                            if (valid_non_edge({a, b})) {
                                // P' = aPb
                                Subgraph P_p({a}, P, {b});

#ifndef NDEBUG
                                // assert that P' is a path of length k
                                assert(P_p.size() == k);
                                for (unsigned i = 0; i < k; ++i)
                                    for (unsigned j = i + 1; j < k; ++j)
                                        if (j - i == 1) assert(valid_edge({P_p[i], P_p[j]}));
                                        else assert(valid_non_edge({P_p[i], P_p[j]}));
#endif

                                if (callback(std::move(P_p), min_vertex)) return true;

                            } else if (valid_edge({a , b}) && with_cycles) { // See 3.4 Listing C_k
                                if (P[0] == min_vertex && P[1] < a) {
                                    // P' = Pba
                                    Subgraph P_p({}, P, {b, a});

#ifndef NDEBUG
                                    // assert that P' is a cycle of length k
                                    assert(P_p.size() == k);
                                    for (unsigned i = 0; i < k; ++i)
                                        for (unsigned j = i + 1; j < k; ++j)
                                            if (j - i == 1 || (i == 0 && j == k-1)) assert(valid_edge({P_p[i], P_p[j]}));
                                            else assert(valid_non_edge({P_p[i], P_p[j]}));
#endif

                                    if (callback(std::move(P_p), min_vertex)) return true;
                                }
                            }
                        }
                    }
                    // no early exit
                    return false;
                }, neighbors, non_neighbors, valid_edge, valid_non_edge);
            } else {
                assert(false);
            }
        }
    };

    template <>
    class CenterFindImpl<3, true> {
    public:
        /**
         *  Find paths and cycles of size 3. Each path or cycle is listed exactly once.
         *
         *  The callback is called with (Subgraph&&, Vertex).
         *
         * @param graph
         * @return
         */
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F neighbors, G /*non_neighbors*/, H valid_edge, I valid_non_edge) {
            Graph::AdjRow Y(graph.size()), Z(graph.size());

            /** P_3, C_3: <y, x, z> **/
            for (Vertex x : graph.vertices()) {
                Y = Z = neighbors(x);
                for (Vertex y : Graph::iterate(Y)) {
                    for (Vertex z : Graph::iterate(Z)) {
                        if (y >= z) continue;
                        assert(y != z); assert(y != x); assert(z != x);

                        if (valid_edge({y, z}) && y < x && y < z && x < z) {
                            // C_3
                            assert(valid_edge({y, x})); assert(valid_edge({y, z})); assert(valid_edge({x, z}));
                            if (callback(Subgraph{y, x, z}, y)) return true;
                        } else if (valid_non_edge({y, z}) && y < z) {
                            // P_3
                            assert(valid_edge({y, x})); assert(valid_non_edge({y, z})); assert(valid_edge({x, z}));
                            Vertex min_vertex = y < x ? y : x;
                            if (callback(Subgraph{y, x, z}, min_vertex)) return true;
                        }
                    }
                }
            }
            return false;
        }
    };

    template <>
    class CenterFindImpl<3, false> {
    public:
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
        /**
         *  Find paths of size 3. Each path is listed exactly once.
         *
         *  The callback is called with (Subgraph&&, Vertex).
         *
         * @param graph
         * @return
         */
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F neighbors, G /*non_neighbors*/, H valid_edge, I valid_non_edge) {
#pragma GCC diagnostic pop
            Graph::AdjRow W(graph.size());

            /** P_3, C_3: <y, x, z> **/
            for (Vertex x : graph.vertices()) {
                W = neighbors(x);
                for (Vertex y : Graph::iterate(W)) {
                    for (Vertex z : Graph::iterate(W)) {
                        if (y >= z) continue;
                        assert(y != z); assert(y != x); assert(z != x);

                        if (valid_non_edge({y, z}) && y < z) {
                            // P_3
                            assert(valid_edge({y, x})); assert(valid_non_edge({y, z})); assert(valid_edge({x, z}));
                            Vertex min_vertex = y < x ? y : x;
                            if (callback(Subgraph{y, x, z}, min_vertex)) return true;
                        }
                    }
                }
            }
            return false;
        }
    };

    template <>
    class CenterFindImpl<2, false> {
    public:
        /**
         * Find paths of size 2. Each path is listed exactly once.
         *
         * The callback is called with (Subgraph&&, Vertex).
         *
         * @param graph
         * @return
         */
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F neighbors, G /*non_neighbors*/, H valid_edge, I /*valid_non_edge*/) {
            /** P_2: <u, v> **/
            for (Vertex u : graph.vertices()) {
                auto V = neighbors(u);
                for (Vertex v : Graph::iterate(V)) {
                    if (u >= v) continue;
                    assert(u != v);

                    if (valid_edge({u, v})) {
                        if (callback(Subgraph{u, v}, u)) return true;
                    }
                }
            }
            return false;
        }
    };

    template <>
    class CenterFindImpl<1, false> {
    public:
        /**
         * Find paths of size 1. Each path is listed exactly once.
         *
         * The callback is called with (Subgraph&&, Vertex).
         *
         * @param graph
         * @return
         */
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F /*neighbors*/, G /*non_neighbors*/, H /*valid_edge*/, I /*valid_non_edge*/) {
            /** P_1: <u> **/
            for (Vertex u : graph.vertices()) {
                if (callback(Subgraph{u}, u)) return true;
            }
            return false;
        }
    };
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_CENTERFINDIMPL_H
