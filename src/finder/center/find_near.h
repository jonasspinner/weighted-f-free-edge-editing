//
// Created by jonas on 31.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FIND_NEAR_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FIND_NEAR_H


#include "find.h"


namespace detail {

    template <int k, bool with_cycles>
    class FindNearImpl {
        static_assert(k >= 4);
    public:
        /**
         * Find paths (and cycles) with length of at least 4 containing the vertices u and v. Each path is listed exactly once.
         *
         * The current implementation is inefficient. It lists all paths (and cycles) with length of at least 4 and filters them.
         *
         * @param graph
         * @param uv
         * @return
         */
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find_near(const Graph& graph, VertexPair uv, SubgraphCallback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            return FindImpl<k, with_cycles>::find(graph, [&](Subgraph &&subgraph) {
                auto vertices = subgraph.vertices();

                bool has_u = std::any_of(vertices.begin(), vertices.end(), [&](Vertex x) { return x == uv.u; });
                bool has_v = std::any_of(vertices.begin(), vertices.end(), [&](Vertex x) { return x == uv.v; });

                if (has_u && has_v) {
                    return callback(std::move(subgraph));
                } else {
                    return false;
                }
            }, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }
    };

    template <>
    class FindNearImpl<3, false> {
    public:
        /**
         * Find paths of length 3 containing the vertices u and v. Each path is listed exactly once.
         *
         * @param graph
         * @param uv
         * @return
         */
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
    };

    template <>
    class FindNearImpl<2, false> {
    public:
        /**
         * Find paths of length 2 containing vertices u and v. Each path is listed exactly once.
         *
         * @param uv
         * @return
         */
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find_near(const Graph& /*graph*/, VertexPair uv, SubgraphCallback callback, F /*neighbors*/, G /*non_neighbors*/, H valid_edge, I /*valid_non_edge*/) {
            assert(false);
            auto [u, v] = uv;

            if (valid_edge(uv)) {
                if (callback(Subgraph{u, v})) return true;
            }
            return false;
        }
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FIND_NEAR_H
