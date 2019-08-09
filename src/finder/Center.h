//
// Created by jonas on 31.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H


#include "../graph/Subgraph.h"
#include "../interfaces/FinderI.h"

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
    class CenterFinderImpl {
    public:
        static_assert(k >= 4);

        // TODO: not 100% sure on implementation details

        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            if constexpr (k >= 6) {

                /** Recursion on finding P_{k-4} **/
                return CenterFinderImpl<k - 4, false>::find(graph, [&](Subgraph&& P) {
                    // P = p_1, ..., p_{k-4}
                    assert(P.size() == k - 4);

                    auto A = neighbors(P[0]); // adjacent to p_1 and non adjacent to P - {p_1}
                    for (int i = 1; i < k - 4; ++i) { A &= non_neighbors(P[i]); }
                    auto B = neighbors(P[k - 5]); // adjacent to p_{k-4} and non adjacent to P - {p_{k-4}}
                    for (int i = 0; i < k - 5; ++i) { B &= non_neighbors(P[i]); }
                    auto C = graph.all_vertices(); // non adjacent to P
                    for (int i = 0; i < k - 4; ++i) { C &= non_neighbors(P[i]); }
#ifndef NDEBUG
                    for (Vertex a : Graph::vertices(A)) {
                        assert(valid_edge({a, P[0]}));
                        for (int i = 1; i < k-4; ++i) assert(valid_non_edge({a, P[i]}));
                    }
                    for (Vertex b : Graph::vertices(B)) {
                        assert(valid_edge({b, P[k-5]}));
                        for (int i = 0; i < k-5; ++i) assert(valid_non_edge({b, P[i]}));
                    }
                    for (Vertex c : Graph::vertices(C)) {
                        for (Vertex p : P) assert(valid_non_edge({c, p}));
                    }
#endif
                    // for each (a, b) \in AxB
                    return Graph::iterate(A, B, [&](Vertex a, Vertex b) {
                        assert(a != b);
                        for (Vertex p : P) { assert(a != p); assert(b != p); }

                        if (valid_non_edge({a, b})) {

                            auto A_p = C & neighbors(a) & non_neighbors(b); // subset of C adjacent to a but not b
                            auto B_p = C & neighbors(b) & non_neighbors(a); // subset of C adjacent to b but not a

                            // for each (u, v) \in A'xB'
                            Graph::iterate(A_p, B_p, [&](Vertex u, Vertex v) {
                                assert(u != v); assert(u != a); assert(u != b); assert(v != a); assert(v != b);
                                for (Vertex p : P) { assert(u != p); assert(v != p); }

                                if (valid_non_edge({u, v})) {
                                    // P' = uaPbv
                                    Subgraph P_p{u, a};
                                    P_p.insert(P_p.end(), P.begin(), P.end());
                                    P_p.push_back(b); P_p.push_back(v);
#ifndef NDEBUG
                                    assert(P_p.size() == k);
                                    for (int i = 0; i < k; ++i)
                                        for (int j = i + 1; j < k; ++j)
                                            if (j - i == 1) assert(valid_edge({P_p[i], P_p[j]}));
                                            else assert(valid_non_edge({P_p[i], P_p[j]}));
#endif
                                    return callback(std::move(P_p));

                                } else if (valid_edge({u, v}) && with_cycles) { // See 3.4 Listing C_k
                                    // P' = Pbvua
                                    Subgraph P_p(P);
                                    P_p.push_back(b); P_p.push_back(v); P_p.push_back(u); P_p.push_back(a);

                                    Vertex min_vertex = P_p[0];
                                    for (int i = 1; i < k; ++i) { min_vertex = std::min(min_vertex, P_p[i]); }

                                    if (P_p[0] == min_vertex && P_p[1] < P_p[k-1]) {
#ifndef NDEBUG
                                        assert(P_p.size() == k);
                                        for (int i = 0; i < k; ++i)
                                            for (int j = i + 1; j < k; ++j)
                                                if (j - i == 1 || j - i == k-1) assert(valid_edge({P_p[i], P_p[j]}));
                                                else assert(valid_non_edge({P_p[i], P_p[j]}));
#endif
                                        return callback(std::move(P_p));
                                    }
                                }
                                return false;
                            });
                        }
                        return false;
                    });
                }, neighbors, non_neighbors, valid_edge, valid_non_edge);
            } else if constexpr (k >= 4) {

                /** Recursion on finding P_{k-2} **/
                return CenterFinderImpl<k - 2, false>::find(graph, [&](Subgraph &&P) {
                    // P = p_1, ..., p_{k-2}
                    assert(P.size() == k - 2);

                    auto A = neighbors(P[0]); // adjacent to p_1 and non adjacent to P - {p_1}
                    for (int i = 1; i < k - 2; ++i) { A &= non_neighbors(P[i]); }
                    auto B = neighbors(P[k - 3]); // adjacent to p_{k-2} and non adjacent to P - {p_{k-2}}
                    for (int i = 0; i < k - 3; ++i) { B &= non_neighbors(P[i]); }
#ifndef NDEBUG
                    for (Vertex a : Graph::vertices(A)) {
                        assert(valid_edge({a, P[0]}));
                        for (int i = 1; i < k-2; ++i) assert(valid_non_edge({a, P[i]}));
                    }
                    for (Vertex b : Graph::vertices(B)) {
                        assert(valid_edge({b, P[k-3]}));
                        for (int i = 0; i < k-3; ++i) assert(valid_non_edge({b, P[i]}));
                    }
#endif
                    // for each (a, b) \in AxB
                    return Graph::iterate(A, B, [&](Vertex a, Vertex b) {
                        assert(a != b);
                        for (Vertex p : P) { assert(a != p); assert(b != p); }

                        if (valid_non_edge({a, b})) {
                            // P' = aPb
                            Subgraph P_p{a};
                            P_p.insert(P_p.end(), P.begin(), P.end());
                            P_p.push_back(b);
#ifndef NDEBUG
                            assert(P_p.size() == k);
                            for (int i = 0; i < k; ++i)
                                for (int j = i + 1; j < k; ++j)
                                    if (j - i == 1) assert(valid_edge({P_p[i], P_p[j]}));
                                    else assert(valid_non_edge({P_p[i], P_p[j]}));
#endif
                            return callback(std::move(P_p));

                        } else if (valid_edge({a , b}) && with_cycles) { // See 3.4 Listing C_k
                            // P' = Pba
                            Subgraph P_p(P);
                            P_p.push_back(b); P_p.push_back(a);

                            Vertex min_vertex = P_p[0];
                            for (int i = 1; i < k; ++i) { min_vertex = std::min(min_vertex, P_p[i]); }

                            if (P_p[0] == min_vertex && P_p[1] < P_p[k-1]) {
#ifndef NDEBUG
                                assert(P_p.size() == k);
                                for (int i = 0; i < k; ++i)
                                    for (int j = i + 1; j < k; ++j)
                                        if (j - i == 1 || (i == 0 && j == k-1)) assert(valid_edge({P_p[i], P_p[j]}));
                                        else assert(valid_non_edge({P_p[i], P_p[j]}));
#endif
                                return callback(std::move(P_p));
                            }
                        }
                        return false;
                    });
                }, neighbors, non_neighbors, valid_edge, valid_non_edge);
            } else {
                assert(false);
            }
        }
    };


    template <>
    class CenterFinderImpl<3, false> {
    public:
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {

            /** P_3: <x, y, z> **/
            for (Vertex x : graph.vertices()) {
                // V - N(x) - {x}
                auto y_candidates = non_neighbors(x);

                for (Vertex y : Graph::vertices(y_candidates)) {
                    for (Vertex z : graph.neighbors(y)) {
                        assert(y != z); assert(y != x); assert(z != x);

                        if (graph.has_edge({x, z}) && y < x) {
                            assert(valid_edge({y, z})); assert(valid_non_edge({y, x})); assert(valid_edge({z, x}));
                            if (callback(Subgraph{y, z, x})) return true;
                        }
                    }
                }
            }
            return false;
        }
    };

    template <>
    class CenterFinderImpl<2, false> {
    public:
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            /** P_2: <u, v> **/
            for (Vertex u : graph.vertices()) {
                for (Vertex v : graph.neighbors(u)) {
                    if (u >= v) continue;
                    assert(u != v);

                    if (valid_edge({u, v})) {
                        if (callback(Subgraph{u, v})) return true;
                    }
                }
            }
            return false;
        }
    };

    template <>
    class CenterFinderImpl<1, false> {
    public:
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            /** P_1: <u> **/
            for (Vertex u : graph.vertices()) {
                if (callback(Subgraph{u})) return true;
            }
            return false;
        }
    };

    template <int length>
    class Center : public FinderI {
        static_assert(length > 1);

    public:
        explicit Center(const Graph &graph) : FinderI(graph) {}

        bool find(SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return  graph.adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair uv) { return  graph.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv); };

            return detail::CenterFinderImpl<length, (length > 3)>::find(graph, callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find(const Graph& forbidden, SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return  graph.adj[u] & ~forbidden.adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.adj[u] & ~forbidden.adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair uv) { return  graph.has_edge(uv) && !forbidden.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv) && !forbidden.has_edge(uv); };

            return detail::CenterFinderImpl<length, (length > 3)>::find(graph, callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find_near(VertexPair uv, SubgraphCallback callback) override { assert(false); return false; }

        bool find_near(VertexPair uv, const Graph& forbidden, SubgraphCallback callback) override  { assert(false); return false; }

    };
}

using CenterRecC4P4 = detail::Center<4>;

#endif //WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H
