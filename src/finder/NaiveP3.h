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
        explicit NaiveP3(const Graph &graph) : FinderI(graph) {}

        bool find(SubgraphCallback callback) override {
            return graph.for_all_vertices([&](Vertex u) {
                return graph.for_all_vertices([&](Vertex v) {
                    if (u >= v) return false;
                    return graph.for_all_vertices([&](Vertex w) {
                        if (v >= w) return false;
                        size_t n_edges = graph.has_edge({u, v}) + graph.has_edge({v, w}) + graph.has_edge({w, u});
                        if (n_edges == 2) {
                            return callback(Subgraph{u, v, w});
                        }
                        return false;
                    });
                });
            });
        }

        bool find(const Graph &forbidden, SubgraphCallback callback) override {
            assert(false);
            return false;
        }

        bool find_near(VertexPair uv, SubgraphCallback callback) override {
            Vertex u = uv.u;
            Vertex v = uv.v;

            return graph.for_all_vertices([&](Vertex w) {
                if (w == u || w == v) return false;
                size_t n_edges = graph.has_edge({u, v}) + graph.has_edge({v, w}) + graph.has_edge({w, u});
                if (n_edges == 2) return callback(Subgraph{u, v, w});
                return false;
            });
        }

        bool find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) override {
            assert(false);
            return false;
        }

    };
}

#endif //CONCEPT_NAIVEP3_H
