//
// Created by jonas on 06.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINTFIND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINTFIND_H


namespace detail {
    template <int k, bool with_cycles, bool intermediate = false>
    class EndpointFindImpl {
        static_assert(k > 2);
    public:
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            Graph::AdjRow A(graph.size()), B(graph.size()), C(graph.size());
            return EndpointFindImpl<k - 2, false, true>::find(graph, [&](Subgraph&& P, Vertex P_min_vertex) {
                assert(P.size() == k - 2);

                C.set();
                for (size_t i = 1; i < k - 3; ++i)
                    C &= non_neighbors(P[i]);

                A = neighbors(P[k-3]); A &= C; A &= non_neighbors(P[0]);
                for (Vertex a : Graph::iterate(A)) {
                    B = neighbors(a); B &= C; B &= non_neighbors(P[k - 3]);

                    for (Vertex b : Graph::iterate(B)) {
                        Vertex min_vertex = std::min({P_min_vertex, a, b});

                        if (valid_non_edge({P[0], b})) {

                            if (P[0] < b || intermediate) {
                                // P' = Pab
                                Subgraph P_p({}, P, {a, b});
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
                            }
                        } else if (valid_edge({P[0], b}) && with_cycles) {

                            if ((P[0] == min_vertex && P[1] < b) || intermediate) {
                                // P' = Pab
                                Subgraph P_p({}, P, {a, b});

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
                return false;
            }, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }
    };

    template <>
    class EndpointFindImpl<3, false, false> {
    public:
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            Graph::AdjRow V(graph.size()), W(graph.size());
            for (Vertex u : graph.vertices()) {
                V = neighbors(u);
                for (Vertex v : Graph::iterate(V)) {
                    W = neighbors(v) & non_neighbors(u);
                    for (Vertex w : Graph::iterate(W))
                        if (u < w) {
                            assert(valid_edge({u, v})); assert(valid_edge({v, w})); assert(valid_non_edge({w, u}));
                            assert(u < w);
                            if (callback(Subgraph{u, v, w}, std::min(u, v)))
                                return true;
                        }
                }
            }
            return false;
        }
    };

    template <>
    class EndpointFindImpl<3, false, true> {
    public:
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F neighbors, G non_neighbors, H /*valid_edge*/, I /*valid_non_edge*/) {
            Graph::AdjRow V(graph.size()), W(graph.size());
            for (Vertex u : graph.vertices()) {
                V = neighbors(u);
                for (Vertex v : Graph::iterate(V)) {
                    W = neighbors(v) & non_neighbors(u);
                    for (Vertex w : Graph::iterate(W))
                        if (callback(Subgraph{u, v, w}, std::min({u, v, w})))
                            return true;
                }
            }
            return false;
        }
    };

    template <>
    class EndpointFindImpl<2, false, false> {
    public:
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F /*neighbors*/, G /*non_neighbors*/, H valid_edge, I /*valid_non_edge*/) {
            for (auto [u, v] : graph.edges()) {
                if (valid_edge({{u, v}}))
                    if (callback(Subgraph{u, v}, u))
                        return true;
            }
            return false;
        }
    };

    template <>
    class EndpointFindImpl<2, false, true> {
    public:
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F neighbors, G /*non_neighbors*/, H /*valid_edge*/, I /*valid_non_edge*/) {
            Graph::AdjRow V(graph.size());
            for (Vertex u : graph.vertices()) {
                V = neighbors(u);
                for (Vertex v : Graph::iterate(V))
                    if (callback(Subgraph{u, v}, std::min(u, v)))
                        return true;
            }
            return false;
        }
    };

    template <>
    class EndpointFindImpl<1, false, true> {
    public:
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F /*neighbors*/, G /*non_neighbors*/, H /*valid_edge*/, I /*valid_non_edge*/) {
            for (Vertex u : graph.vertices()) {
                if (callback(Subgraph{u}, u))
                    return true;
            }
            return false;
        }
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINTFIND_H
