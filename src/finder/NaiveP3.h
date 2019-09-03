//
// Created by jonas on 02.07.19.
//

#ifndef CONCEPT_NAIVEP3_H
#define CONCEPT_NAIVEP3_H


#include "../graph/Graph.h"
#include "../graph/Subgraph.h"

namespace Finder {
    class NaiveP3 : public FinderI {

    public:
        explicit NaiveP3(const Graph &graph_ref) : FinderI(graph_ref) {}

        /**
         * Calls callback for all P_3's.
         *
         * @param callback
         * @return
         */
        bool find(SubgraphCallback callback) override {
            for (Vertex u : graph.vertices()) {
                for (Vertex v : graph.vertices()) {
                    if (u >= v) continue;
                    for (Vertex w : graph.vertices()) {
                        if (v >= w) continue;
                        int n_edges = graph.has_edge({u, v}) + graph.has_edge({v, w}) + graph.has_edge({w, u});
                        if (n_edges == 2) {
                            if (callback(Subgraph{u, v, w})) return true;
                        }
                    }
                }
            }
            return false;
        }

        /**
         * Calls callback for all P_3's. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
         *
         * @param forbidden
         * @param callback
         * @return
         */
        bool find(const Graph &forbidden, SubgraphCallback callback) override {
            return find([&](Subgraph &&subgraph) {
                Vertex u = subgraph[0], v = subgraph[1], w = subgraph[2];
                if (forbidden.has_edge({u, v}) || forbidden.has_edge({u, w}) || forbidden.has_edge({v, w})) return false;
                return callback(std::move(subgraph));
            });
        }

        /**
         * Calls callback for all P_3's having both u and v as vertices.
         *
         * @param uv
         * @param callback
         * @return
         */
        bool find_near(VertexPair uv, SubgraphCallback callback) override {
            auto [u, v] = uv;

            for (Vertex w : graph.vertices()) {
                if (w == u || w == v) continue;
                int n_edges = graph.has_edge({u, v}) + graph.has_edge({v, w}) + graph.has_edge({w, u});
                if (n_edges == 2) {
                    if (callback(Subgraph{u, v, w})) return true;
                }
            }
            return false;
        }

        /**
         * Calls callback for all P_3's having both u and v as vertices. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
         *
         * @param uv
         * @param forbidden
         * @param callback
         * @return
         */
        bool find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) override {
            return find_near(uv, [&](Subgraph &&subgraph) {
                Vertex u = subgraph[0], v = subgraph[1], w = subgraph[2];
                if (forbidden.has_edge({u, v}) || forbidden.has_edge({u, w}) || forbidden.has_edge({v, w})) return false;
                return callback(std::move(subgraph));
            });
        }

    };
}

#endif //CONCEPT_NAIVEP3_H
