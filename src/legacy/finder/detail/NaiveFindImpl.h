//
// Created by jonas on 02.11.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_NAIVEFINDIMPL_H
#define WEIGHTED_F_FREE_EDGE_EDITING_NAIVEFINDIMPL_H


namespace detail {
    template <int k>
    class OrderedVertexSetIterator {
    public:
        template <typename VertexSetCallback>
        static bool iterate(const Graph& graph, VertexSetCallback callback) {
            return OrderedVertexSetIterator<k - 1>::iterate(graph, [&](Subgraph&& P) {
                for (Vertex a : graph.vertices()) {
                    if (std::all_of(P.vertices().begin(), P.vertices().end(), [&](auto u) { return u != a; })) {
                        Subgraph P_p({}, P, {a});
                        if (callback(std::move(P_p))) return true;
                    }
                }
                return false;
            });
        }
    };

    template <>
    class OrderedVertexSetIterator<1> {
    public:
        template <typename VertexSetCallback>
        static bool iterate(const Graph& graph, VertexSetCallback callback) {
            for (Vertex a : graph.vertices()) {
                if (callback({a})) return true;
            }
            return false;
        }
    };

    template <int k, bool with_cycles>
    class NaiveFindImpl {
        static_assert(k > 2);
    public:
        template <typename SubgraphCallback, typename H, typename I>
        static bool find(const Graph& graph, SubgraphCallback callback, H valid_edge, I valid_non_edge) {
            Graph::AdjRow A(graph.size()), B(graph.size()), C(graph.size());
            return OrderedVertexSetIterator<k>::iterate(graph, [&](Subgraph&& P) {
                assert(P.size() == k);

                if (valid_non_edge({P[0], P[k-1]})) {
                    if (P[0] < P[k-1]) {
                        // check if P' is a path of length k
                        bool is_pk = true;
                        for (unsigned i = 0; i < k; ++i)
                            for (unsigned j = i + 1; j < k; ++j)
                                if (j - i == 1) {
                                    is_pk &= valid_edge({P[i], P[j]});
                                } else {
                                    is_pk &= valid_non_edge({P[i], P[j]});
                                }

                        if (is_pk)
                            if (callback(std::move(P))) return true;
                    }
                } else if (valid_edge({P[0], P[k-1]}) && with_cycles) {
                    Vertex min_vertex = *std::min_element(P.vertices().begin(), P.vertices().end());
                    if ((P[0] == min_vertex && P[1] < P[k-1])) {

                        // check if P' is a cycle of length k
                        bool is_ck = true;
                        for (unsigned i = 0; i < k; ++i)
                            for (unsigned j = i + 1; j < k; ++j)
                                if (j - i == 1 || j - i == k-1) {
                                    is_ck &= valid_edge({P[i], P[j]});
                                } else {
                                    is_ck &= valid_non_edge({P[i], P[j]});
                                }

                        if (is_ck)
                            if (callback(std::move(P))) return true;
                    }
                }
                return false;
            });
        }
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_NAIVEFINDIMPL_H
