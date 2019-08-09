//
// Created by jonas on 23.07.19.
//

#ifndef CONCEPT_CENTERP3_H
#define CONCEPT_CENTERP3_H


#include "../interfaces/FinderI.h"

namespace Finder {
    class CenterP3 : public FinderI {

        Graph::AdjRow w_candidates;

    public:
        explicit CenterP3(const Graph &graph) : FinderI(graph), w_candidates(graph.size()) {}

        bool find(SubgraphCallback callback) override {
            auto valid_edge = [&](VertexPair uv) { return graph.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv); };

            return graph.for_all_vertices([&](Vertex u) {
                auto v_candidates = ~graph.adj[u]; v_candidates[u] = false;

                return Graph::iterate(v_candidates, [&](Vertex v) {
                    return graph.for_neighbors_of(v, [&](Vertex w) {
                        if (graph.has_edge({u, w}) && v < u) {

                            assert(valid_edge({u, w})); assert(valid_edge({v, w})); assert(valid_non_edge({u, v}));
                            return callback(Subgraph{u, v, w});
                        } else { return false; }
                    });
                });
            });
        }

        bool find(const Graph &forbidden, SubgraphCallback callback) override {
            auto neighbors = [&](Vertex u) { return graph.adj[u] & ~forbidden.adj[u]; };
            auto non_neighbors = [&](Vertex u) { return ~graph.adj[u] & ~forbidden.adj[u]; };
            auto valid_edge = [&](VertexPair uv) { return graph.has_edge(uv) && !forbidden.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv) && !forbidden.has_edge(uv); };

            return graph.for_all_edges([&](VertexPair uv) {
                if (forbidden.has_edge(uv)) return false;
                Vertex u = uv.u;
                Vertex v = uv.v;

                w_candidates = neighbors(u) & non_neighbors(v);
                w_candidates[v] = false;

                bool exited = Graph::iterate(w_candidates, [&](Vertex w) {
                    assert(valid_edge({u, v})); assert(valid_edge({u, w})); assert(valid_non_edge({v, w}));
                    return callback(Subgraph{u, v, w});
                });
                if (exited) return true;

                w_candidates = neighbors(v) & non_neighbors(u);
                w_candidates[u] = false;

                exited = Graph::iterate(w_candidates, [&](Vertex w) {
                    assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_edge({v, w}));
                    return callback(Subgraph{v, u, w});
                });

                return exited;
            });
        }

        bool find_near(VertexPair uv, SubgraphCallback callback) override {
            Vertex u = uv.u;
            Vertex v = uv.v;

            if (graph.has_edge(uv)) {
                w_candidates = graph.adj[u] & ~graph.adj[v];
                w_candidates[v] = false;

                bool exited = Graph::iterate(w_candidates, [&](Vertex w) {
                    return callback(Subgraph{u, v, w});
                });
                if (exited) return true;

                w_candidates = ~graph.adj[u] & graph.adj[v];
                w_candidates[u] = false;

                exited = Graph::iterate(w_candidates, [&](Vertex w) {
                    return callback(Subgraph{v, u, w});
                });

                return exited;
            } else {
                bool exited = Graph::iterate(graph.adj[u] & graph.adj[v], [&](Vertex w) {
                    return callback(Subgraph{w, u, v});
                });

                return exited;
            }
        }

        bool find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) override {
            auto neighbors = [&](Vertex u) { return graph.adj[u] & ~forbidden.adj[u]; };
            auto non_neighbors = [&](Vertex u) { return ~graph.adj[u] & ~forbidden.adj[u]; };
            auto valid_edge = [&](VertexPair uv) { return graph.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv); };
            Vertex u = uv.u;
            Vertex v = uv.v;

            if (valid_edge(uv)) {
                w_candidates = neighbors(u) & non_neighbors(v);
                w_candidates[v] = false;

                bool exited = Graph::iterate(w_candidates, [&](Vertex w) {
                    return callback(Subgraph{u, v, w});
                });
                if (exited) return true;

                w_candidates = neighbors(v) & non_neighbors(u);
                w_candidates[u] = false;

                exited = Graph::iterate(w_candidates, [&](Vertex w) {
                    return callback(Subgraph{v, u, w});
                });

                return exited;
            } else if (valid_non_edge(uv)) {

                w_candidates = neighbors(u) & neighbors(v);
                bool exited = Graph::iterate(w_candidates, [&](Vertex w) {
                    return callback(Subgraph{w, u, v});
                });

                return exited;
            }
            return false;
        }

    };
}


#endif //CONCEPT_CENTERP3_H
