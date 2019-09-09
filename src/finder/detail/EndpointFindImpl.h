//
// Created by jonas on 06.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINTFIND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINTFIND_H


namespace detail {
    template <int k, bool with_cycles>
    class EndpointFindImpl {
        static_assert(k > 2);
    public:
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            Graph::AdjRow A(graph.size()), B(graph.size()), C(graph.size());
            return EndpointFindImpl<k - 2, with_cycles>::find(graph, [&](Subgraph&& P, Vertex P_min_vertex) {
                assert(P.size() == k - 2);

                C = Graph::AdjRow(graph.size());
                C.set();
                for (size_t i = 1; i < k - 3; ++i) {
                    C -= neighbors(P[i]);
                    C[P[i]] = false;
                }

                A = neighbors(P[k-3]); A &= C; A -= neighbors(P[0]); A[P[0]] = false;
                for (Vertex a : Graph::iterate(A)) {
                    B = neighbors(a); B &= C; B -= neighbors(P[k-3]); B[P[k-3]] = false;

                    for (Vertex b : Graph::iterate(B)) {
                        if (valid_non_edge({P[0], b}) && P[0] < b) {
                            // P' = Pab
                            Subgraph P_p(P);
                            P_p.push_back(a); P_p.push_back(b);

                            Vertex min_vertex = std::min({P_min_vertex, a, b});
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
                        } else if (valid_edge({P[0], b}) && with_cycles) {
                            // P' = Pab
                            Subgraph P_p(P);
                            P_p.push_back(a); P_p.push_back(b);

                            Vertex min_vertex = std::min({P_min_vertex, a, b});

                            if (P_p[0] == min_vertex && P[1] < P[k]) {

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
    class EndpointFindImpl<2, false> {
    public:
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, F /*neighbors*/, G /*non_neighbors*/, H valid_edge, I /*valid_non_edge*/) {
            for (VertexPair uv : graph.edges()) {
                if (valid_edge(uv))
                    if (callback(Subgraph{uv.u, uv.v}, uv.u))
                        return true;
            }
            return false;
        }
    };

    template <>
    class EndpointFindImpl<1, false> {
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
