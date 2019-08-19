//
// Created by jonas on 31.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H


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
            Graph::AdjRow A(graph.size()), B(graph.size());

            if constexpr (k >= 6) {
                Graph::AdjRow C(graph.size()), A_p(graph.size()), B_p(graph.size());

                /** Recursion on finding P_{k-4} **/
                return CenterFinderImpl<k - 4, false>::find(graph, [&](Subgraph&& P) {
                    // P = p_1, ..., p_{k-4}
                    assert(P.size() == k - 4);

                    A = neighbors(P[0]); // adjacent to p_1 and non adjacent to P - {p_1}
                    for (unsigned i = 1; i < k - 4; ++i) { A &= non_neighbors(P[i]); }

                    B = neighbors(P[k - 5]); // adjacent to p_{k-4} and non adjacent to P - {p_{k-4}}
                    for (unsigned i = 0; i < k - 5; ++i) { B &= non_neighbors(P[i]); }

                    C = graph.all_vertices(); // non adjacent to P
                    for (unsigned i = 0; i < k - 4; ++i) { C &= non_neighbors(P[i]); }

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

                            assert(a != b);
                            for (Vertex p : P.vertices()) { assert(a != p); assert(b != p); }

                            if (valid_non_edge({a, b})) {

                                A_p = C & neighbors(a) & non_neighbors(b); // subset of C adjacent to a but not b
                                B_p = C & neighbors(b) & non_neighbors(a); // subset of C adjacent to b but not a

                                // for each (u, v) \in A'xB'
                                for (Vertex u : Graph::iterate(A_p)) {
                                    for (Vertex v : Graph::iterate(B_p)) {

                                        assert(u != v); assert(u != a); assert(u != b); assert(v != a); assert(v != b);
                                        for (Vertex p : P.vertices()) { assert(u != p); assert(v != p); }

                                        if (valid_non_edge({u, v})) {
                                            // P' = uaPbv
                                            Subgraph P_p{u, a};
                                            P_p.append(P);
                                            P_p.push_back(b); P_p.push_back(v);
#ifndef NDEBUG
                                            assert(P_p.size() == k);
                                            for (unsigned i = 0; i < k; ++i)
                                                for (unsigned j = i + 1; j < k; ++j)
                                                    if (j - i == 1) {
                                                        assert(valid_edge({P_p[i], P_p[j]}));
                                                    } else {
                                                        assert(valid_non_edge({P_p[i], P_p[j]}));
                                                    }
#endif
                                            if (callback(std::move(P_p))) return true;

                                        } else if (valid_edge({u, v}) && with_cycles) { // See 3.4 Listing C_k
                                            // P' = Pbvua
                                            Subgraph P_p(P);
                                            P_p.push_back(b); P_p.push_back(v); P_p.push_back(u); P_p.push_back(a);

                                            auto vertices = P_p.vertices();
                                            Vertex min_vertex = *std::min_element(vertices.begin(), vertices.end());

                                            if (P_p[0] == min_vertex && P_p[1] < P_p[k-1]) {
#ifndef NDEBUG
                                                assert(P_p.size() == k);
                                                for (unsigned i = 0; i < k; ++i)
                                                    for (unsigned j = i + 1; j < k; ++j)
                                                        if (j - i == 1 || j - i == k-1) {
                                                            assert(valid_edge({P_p[i], P_p[j]}));
                                                        } else {
                                                            assert(valid_non_edge({P_p[i], P_p[j]}));
                                                        }
#endif
                                                if (callback(std::move(P_p))) return true;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    return false;
                }, neighbors, non_neighbors, valid_edge, valid_non_edge);
            } else if constexpr (k >= 4) {

                /** Recursion on finding P_{k-2} **/
                return CenterFinderImpl<k - 2, false>::find(graph, [&](Subgraph &&P) {
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
                    return Graph::iterate(A, B, [&](Vertex a, Vertex b) {
                        assert(a != b);
                        for (Vertex p : P.vertices()) { assert(a != p); assert(b != p); }

                        if (valid_non_edge({a, b})) {
                            // P' = aPb
                            Subgraph P_p{a};
                            P_p.append(P);
                            P_p.push_back(b);
#ifndef NDEBUG
                            assert(P_p.size() == k);
                            for (unsigned i = 0; i < k; ++i)
                                for (unsigned j = i + 1; j < k; ++j)
                                    if (j - i == 1) assert(valid_edge({P_p[i], P_p[j]}));
                                    else assert(valid_non_edge({P_p[i], P_p[j]}));
#endif
                            return callback(std::move(P_p));

                        } else if (valid_edge({a , b}) && with_cycles) { // See 3.4 Listing C_k
                            // P' = Pba
                            Subgraph P_p(P);
                            P_p.push_back(b); P_p.push_back(a);

                            auto vertices = P_p.vertices();
                            Vertex min_vertex = *std::min_element(vertices.begin(), vertices.end());

                            if (P_p[0] == min_vertex && P_p[1] < P_p[k-1]) {
#ifndef NDEBUG
                                assert(P_p.size() == k);
                                for (unsigned i = 0; i < k; ++i)
                                    for (unsigned j = i + 1; j < k; ++j)
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
            Graph::AdjRow Y(graph.size()), Z(graph.size());

            /** P_3: <x, y, z> **/
            for (Vertex x : graph.vertices()) {
                // V - N(x) - {x}
                Y = non_neighbors(x);
                for (Vertex y : Graph::iterate(Y)) {
                    if (y >= x) continue;

                    Z = neighbors(y) & neighbors(x);
                    for (Vertex z : Graph::iterate(Z)) {
                        assert(y != z); assert(y != x); assert(z != x);

                        assert(valid_edge({y, z})); assert(valid_non_edge({y, x})); assert(valid_edge({z, x}));
                        if (callback(Subgraph{y, z, x})) return true;
                    }
                }
            }
            return false;
        }

        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find_near(const Graph& graph, VertexPair uv, SubgraphCallback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            assert(false);
            Graph::AdjRow Z(graph.size());
            auto [u, v] = uv;

            if (valid_edge(uv)) {
                Z = neighbors(v) & non_neighbors(u);
                for (Vertex z : Graph::iterate(Z)) {
                    assert(valid_edge({u, v})); assert(valid_non_edge({u, z})); assert(valid_edge({v, z}));
                    if (callback(Subgraph{u, v, z})) return true;
                }

                Z = neighbors(u) & non_neighbors(v);
                for (Vertex z : Graph::iterate(Z)) {
                    assert(valid_edge({z, u})); assert(valid_non_edge({z, v})); assert(valid_edge({u, v}));
                    if (callback(Subgraph{z, u, v})) return true;
                }
            } else if (valid_non_edge(uv)) {
                Z = neighbors(uv.u) & neighbors(uv.v);
                for (Vertex z : Graph::iterate(Z)) {
                    assert(valid_edge({u, z})); assert(valid_non_edge({u, v})); assert(valid_edge({z, v}));
                    if (callback(Subgraph{u, z, v})) return true;
                }
            }
            return false;
        }

        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find_near(const Graph& graph, Vertex u, SubgraphCallback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            assert(false);
            Graph::AdjRow B(graph.size());

            auto A = neighbors(u);
            for (Vertex a : Graph::iterate(A)) {

                B = neighbors(a) & non_neighbors(u);
                for (Vertex b : Graph::iterate(B)) {
                    if (u < b) {
                        assert(valid_edge({u, a})); assert(valid_non_edge({u, b})); assert(valid_edge({a, b}));
                        if (callback(Subgraph{u, a, b})) return true;
                    } else {
                        assert(valid_edge({b, a})); assert(valid_non_edge({b, u})); assert(valid_edge({a, u}));
                        if (callback(Subgraph{b, a, u})) return true;
                    }
                }

                B = neighbors(u);
                for (Vertex b : Graph::iterate(B)) {
                    if (a == b) {
                        continue;
                    } else if (a < b) {
                        assert(valid_edge({a, u})); assert(valid_non_edge({a, b})); assert(valid_edge({u, b}));
                        if (callback(Subgraph{a, u, b})) return true;
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
        static bool find(const Graph& graph, SubgraphCallback callback, F neighbors, G /*non_neighbors*/, H valid_edge, I /*valid_non_edge*/) {
            /** P_2: <u, v> **/
            for (Vertex u : graph.vertices()) {
                auto V = neighbors(u);
                for (Vertex v : Graph::iterate(V)) {
                    if (u >= v) continue;
                    assert(u != v);

                    if (valid_edge({u, v})) {
                        if (callback(Subgraph{u, v})) return true;
                    }
                }
            }
            return false;
        }

        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find_near(const Graph& /*graph*/, VertexPair uv, SubgraphCallback callback, F /*neighbors*/, G /*non_neighbors*/, H valid_edge, I /*valid_non_edge*/) {
            assert(false);
            auto [u, v] = uv;

            if (valid_edge(uv)) {
                if (callback(Subgraph{u, v})) return true;
            }
            return false;
        }

        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find_near(const Graph& /*graph*/, Vertex u, SubgraphCallback callback, F neighbors, G /*non_neighbors*/, H valid_edge, I /*valid_non_edge*/) {
            assert(false);
            auto Z = neighbors(u);
            for (Vertex z : Graph::iterate(Z)) {
                if (u < z) {
                    assert(valid_edge({u, z}));
                    if (callback(Subgraph{u, z})) return true;
                } else {
                    assert(valid_edge({u, z}));
                    if (callback(Subgraph{z, u})) return true;
                }
            }
            return false;
        }
    };

    template <>
    class CenterFinderImpl<1, false> {
    public:
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F /*neighbors*/, G /*non_neighbors*/, H /*valid_edge*/, I /*valid_non_edge*/) {
            /** P_1: <u> **/
            for (Vertex u : graph.vertices()) {
                if (callback(Subgraph{u})) return true;
            }
            return false;
        }

        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find_near(const Graph& /*graph*/, Vertex u, SubgraphCallback callback, F /*neighbors*/, G /*non_neighbors*/, H /*valid_edge*/, I /*valid_non_edge*/) {
            assert(false);
            return callback(Subgraph{u});
        }
    };

    template <int length>
    class Center : public FinderI {
        static_assert(length > 1);

    public:
        explicit Center(const Graph &graph_ref) : FinderI(graph_ref) {}

        bool find(SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return  graph.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair uv) { return  graph.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv); };

            return detail::CenterFinderImpl<length, (length > 3)>::find(graph, callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find(const Graph& forbidden, SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return graph.m_adj[u] & ~forbidden.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u] & ~forbidden.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair uv) { return  graph.has_edge(uv) && !forbidden.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv) && !forbidden.has_edge(uv); };

            return detail::CenterFinderImpl<length, (length > 3)>::find(graph, callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find_near(VertexPair uv, SubgraphCallback callback) override {
            return find([&](Subgraph &&subgraph) {
                auto vertices = subgraph.vertices();

                bool has_u = std::any_of(vertices.begin(), vertices.end(), [&](Vertex x) { return x == uv.u; });
                bool has_v = std::any_of(vertices.begin(), vertices.end(), [&](Vertex x) { return x == uv.v; });

                if (has_u && has_v) {
                    return callback(std::move(subgraph));
                } else {
                    return false;
                }
            });
        }

        bool find_near(VertexPair /*uv*/, const Graph& /*forbidden*/, SubgraphCallback /*callback*/) override  { assert(false); return false; }

    };
}

using CenterRecC4P4 = detail::Center<4>;

#endif //WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H
