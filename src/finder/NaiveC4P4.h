//
// Created by jonas on 29.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_NAIVEC4P4_H
#define WEIGHTED_F_FREE_EDGE_EDITING_NAIVEC4P4_H


#include "../interfaces/FinderI.h"

namespace Finder {
    class NaiveC4P4 : public FinderI {
    public:
        explicit NaiveC4P4(const Graph& graph_ref) : FinderI(graph_ref) {}

        template <typename H, typename I>
        bool find(const SubgraphCallback& callback, H valid_edge, I valid_non_edge) {

            for (Vertex u : graph.vertices()) {
                for (Vertex v : graph.vertices()) {
                    for (Vertex w : graph.vertices()) {
                        for (Vertex x : graph.vertices()) {
                            if (u != v && u != w && u != x && v != w && v != x && w != x) {
                                /* u-v-w-x */
                                if (valid_edge({u, v}) && valid_non_edge({u, w}) && valid_edge({v, w}) && valid_non_edge({v, x}) && valid_edge({w, x})) {
                                    if (valid_non_edge({u, x})) {
                                        // P_4
                                        if (u < x) // p_1 < p_k
                                            if (callback(Subgraph({u, v, w, x}))) return true;
                                    } else if (valid_edge({u, x})) {
                                        // C_4
                                        if (u < std::min({v, w, x})) // p_1 is smallest
                                            if (v < x) // p_2 < p_k
                                                if (callback(Subgraph({u, v, w, x}))) return true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return false;
        }

        bool find(SubgraphCallback callback) override {
            auto valid_edge =     [&](VertexPair uv) { return graph.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv); };

            return find(callback, valid_edge, valid_non_edge);
        }

        bool find(const Graph &forbidden, SubgraphCallback callback) override {
            auto valid_edge =     [&](VertexPair uv) { return graph.has_edge(uv) && !forbidden.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv) && !forbidden.has_edge(uv); };

            return find(callback, valid_edge, valid_non_edge);
        }

        template <typename H, typename I>
        bool find_near(VertexPair uv, const SubgraphCallback& callback, H valid_edge, I valid_non_edge) {

            return find([&](Subgraph &&subgraph) {
                auto vertices = subgraph.vertices();

                bool has_u = std::any_of(vertices.begin(), vertices.end(), [&](Vertex x) { return x == uv.u; });
                bool has_v = std::any_of(vertices.begin(), vertices.end(), [&](Vertex x) { return x == uv.v; });

                if (has_u && has_v) {
                    return callback(std::move(subgraph));
                } else {
                    return false;
                }
            }, valid_edge, valid_non_edge);
        }

        bool find_near(VertexPair uv, SubgraphCallback callback) override {
            auto valid_edge = [&](VertexPair xy) { return graph.has_edge(xy); };
            auto valid_non_edge = [&](VertexPair xy) { return !graph.has_edge(xy); };

            return find_near(uv, callback, valid_edge, valid_non_edge);
        };

        bool find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) override {
            auto valid_edge = [&](VertexPair xy) { return graph.has_edge(xy) && !forbidden.has_edge(xy); };
            auto valid_non_edge = [&](VertexPair xy) { return !graph.has_edge(xy) && !forbidden.has_edge(xy); };

            return find_near(uv, callback, valid_edge, valid_non_edge);
        }

    };
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_NAIVEC4P4_H
