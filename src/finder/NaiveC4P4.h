//
// Created by jonas on 29.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_NAIVEC4P4_H
#define WEIGHTED_F_FREE_EDGE_EDITING_NAIVEC4P4_H


#include "../interfaces/FinderI.h"

namespace Finder {
    class NaiveC4P4 : public FinderI {
    public:
        explicit NaiveC4P4(const Graph& graph) : FinderI(graph) {}

        bool find(SubgraphCallback callback) override {
            /* w?x
             * | |
             * u-v
             */
            /*
            return graph.for_all_vertices([&](Vertex u) {
                return graph.for_all_vertices([&](Vertex v) {
                    if ((u == v) || !graph.has_edge({u, v})) return false;
                    return graph.for_all_vertices([&](Vertex w) {
                        if ((u == w) || (v == w) || !graph.has_edge({u, w}) || graph.has_edge({v, w})) return false;
                        return graph.for_all_vertices([&](Vertex x) {
                            if ((u == x) || (v == x) || (w == x) || graph.has_edge({u, x}) || !graph.has_edge({v, x})) return false;
                            if (graph.has_edge({w, x}) && (u >= v || u >= w || u >= x)) return false;
                            return callback(Subgraph({u, v, w, x}));
                        });
                    });
                });
            });*/

            /* x-w
             * ? |
             * u-v */
            return graph.for_all_vertices([&](Vertex u) {
                return graph.for_neighbors_of(u, [&](Vertex v) {
                    // if (u >= v) return false;
                    return graph.for_neighbors_of(v, [&](Vertex w) {
                        if (u == w || graph.has_edge({u, w})) return false;
                        return graph.for_neighbors_of(w, [&](Vertex x) {
                            if (v == x || graph.has_edge({v, x})) return false;
                            if (u > x) return false;
                            if (graph.has_edge({w, x}) && (u > w || w > x)) return false;
                            // if (graph.has_edge({w, x}) && (u >= v || u >= w || u >= x)) return false;
                            return callback(Subgraph({u, v, w, x}));
                        });
                    });
                });
            });
        }

        bool find(const Graph &forbidden, SubgraphCallback callback) override {
            return find([&](Subgraph&& subgraph) {
                bool touched = subgraph.for_all_vertex_pairs([&](VertexPair uv) {
                    return forbidden.has_edge(uv);
                });
                return !touched && callback(std::move(subgraph));
            });
        }

        bool find_near(VertexPair uv, SubgraphCallback callback) override { assert(false); return false; };

        bool find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) override { assert(false); return false; }

    };
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_NAIVEC4P4_H
