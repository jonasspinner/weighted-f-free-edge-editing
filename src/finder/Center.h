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
        // static_assert(!with_cycles);

        template <typename SubgraphCallback>
        static bool find(const Graph& graph, SubgraphCallback callback) {
            if constexpr (k >= 6) {
                /** Recursion on finding P_{k-4} **/
                return CenterFinderImpl<k - 4, false>::find(graph, [&](Subgraph&& P) {
                    // P = p_1, ..., p_{k-4}
                    assert(P.size() == k - 4);

                    auto A = graph.adj[P[0]]; // adjacent to p_1 and non adjacent to P - {p_1}
                    for (int i = 1; i < k - 4; ++i) { A &= ~graph.adj[P[i]]; }
                    auto B = graph.adj[P[k - 5]]; // adjacent to p_{k-4} and non adjacent to P - {p_{k-4}}
                    for (int i = 0; i < k - 5; ++i) { B &= ~graph.adj[P[i]]; }
                    auto C = graph.all_vertices(); // non adjacent to P
                    for (int i = 0; i < k - 4; ++i) { C &= ~graph.adj[P[i]]; }

                    // for each (a, b) \in AxB
                    return Graph::iterate(A, B, [&](Vertex a, Vertex b) {
                        if (!graph.has_edge({a, b})) {

                            auto A_p = C & graph.adj[a] & ~graph.adj[b]; // subset of C adjacent to a but not b
                            auto B_p = C & graph.adj[b] & ~graph.adj[a]; // subset of C adjacent to b but not a

                            // for each (u, v) \in A'xB'
                            Graph::iterate(A_p, B_p, [&](Vertex u, Vertex v) {
                                if (!graph.has_edge({u, v})) {
                                    // P' = uaPbv
                                    Subgraph P_p{u, a};
                                    P_p.insert(P_p.end(), P.begin(), P.end());
                                    P_p.push_back(b); P_p.push_back(v);
                                    return callback(std::move(P_p));
                                } else if constexpr (with_cycles) { // See 3.4 Listing C_k
                                    // P' = Pbvua
                                    Subgraph P_p(P);
                                    P_p.push_back(b); P_p.push_back(v); P_p.push_back(u); P_p.push_back(a);
                                    Vertex min_vertex = P_p[0];
                                    for (int i = 1; i < k; ++i) { min_vertex = std::min(min_vertex, P_p[i]); }
                                    if (P_p[0] == min_vertex && P_p[1] < P_p[k-1]) {
                                        return callback(std::move(P_p));
                                    }
                                }
                                return false;
                            });
                        }
                        return false;
                    });
                });
            } else if constexpr (k >= 4) {
                /** Recursion on finding P_{k-2} **/
                return CenterFinderImpl<k - 2, false>::find(graph, [&](Subgraph &&P) {
                    // P = p_1, ..., p_{k-2}
                    assert(P.size() == k - 2);

                    auto A = graph.adj[P[0]]; // adjacent to p_1 and non adjacent to P - {p_1}
                    for (int i = 1; i < k - 2; ++i) { A &= ~graph.adj[P[i]]; }
                    auto B = graph.adj[P[k - 3]]; // adjacent to p_{k-2} and non adjacent to P - {p_{k-2}}
                    for (int i = 0; i < k - 3; ++i) { B &= ~graph.adj[P[i]]; }

                    // for each (a, b) \in AxB
                    return Graph::iterate(A, B, [&](Vertex a, Vertex b) {
                        if (!graph.has_edge({a, b})) {
                            // P' = aPb
                            Subgraph P_p{a};
                            P_p.insert(P_p.end(), P.begin(), P.end());
                            P_p.push_back(b);
                            return callback(std::move(P_p));
                        } else if constexpr (with_cycles) { // See 3.4 Listing C_k
                            // P' = Pba
                            Subgraph P_p(P);
                            P_p.push_back(b); P_p.push_back(a);
                            Vertex min_vertex = P_p[0];
                            for (int i = 1; i < k; ++i) { min_vertex = std::min(min_vertex, P_p[i]); }
                            if (P_p[0] == min_vertex && P_p[1] < P_p[k-1]) {
                                return callback(std::move(P_p));
                            }
                        }
                        return false;
                    });
                });
            } else {
                assert(false);
            }
        }
    };


    template <>
    class CenterFinderImpl<3, false> {
    public:
        template <typename SubgraphCallback>
        static bool find(const Graph& graph, SubgraphCallback callback) {
            /** P_3: <x, y, z> **/
            return graph.for_all_vertices([&](Vertex x) {
                auto y_candidates = ~graph.adj[x]; // V - N(x) - {x}
                y_candidates[x] = false;
                return Graph::iterate(y_candidates, [&](Vertex y) {
                    return graph.for_neighbors_of(y, [&](Vertex z) {
                        if (graph.has_edge({x, z}) && y < x) {
                            return callback(Subgraph{y, z, x});
                        }
                        return false;
                    });
                });
            });
        }
    };

    template <>
    class CenterFinderImpl<2, false> {
    public:
        template <typename SubgraphCallback>
        static bool find(const Graph& graph, SubgraphCallback callback) {
            /** P_2: <u, v> **/
            return graph.for_all_vertices([&](Vertex u) {
                return graph.for_all_vertices([&](Vertex v) {
                    if (u >= v) return false;
                    return callback(Subgraph{u, v});
                });
            });
        }
    };

    template <>
    class CenterFinderImpl<1, false> {
    public:
        template <typename SubgraphCallback>
        static bool find(const Graph& graph, SubgraphCallback callback) {
            /** P_1: <u> **/
            return graph.for_all_vertices([&](Vertex u) {
                return callback(Subgraph{u});
            });
        }
    };

    template <int length>
    class Center : public FinderI {
        static_assert(length > 1);

    public:
        explicit Center(const Graph &graph) : FinderI(graph) {}

        bool find(SubgraphCallback callback) override {
            return detail::CenterFinderImpl<length, (length > 3)>::find(graph, callback);
        }

        bool find(const Graph& forbidden, SubgraphCallback callback) override { assert(false); return false; }

        bool find_near(VertexPair uv, SubgraphCallback callback) override { assert(false); return false; }

        bool find_near(VertexPair uv, const Graph& forbidden, SubgraphCallback callback) override  { assert(false); return false; }

    };
}

using CenterRecC4P4 = detail::Center<4>;

#endif //WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H
