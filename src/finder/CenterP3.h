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
        explicit CenterP3(const Graph &graph) : FinderI(graph), w_candidates(graph.n_vertices()) {}

        bool find(SubgraphCallback callback) override {

            return graph.for_all_edges([&](VertexPair uv) {
                Vertex u = uv.u;
                Vertex v = uv.v;

                w_candidates = graph.adj[u] & ~graph.adj[v];
                w_candidates[v] = false;

                bool exited = Graph::iterate(w_candidates, [&](Vertex w) {
                    return callback(Subgraph{u, v, w});
                });

                w_candidates = ~graph.adj[u] & graph.adj[v];
                w_candidates[u] = false;

                exited |= Graph::iterate(w_candidates, [&](Vertex w) {
                    return callback(Subgraph{v, u, w});
                });

                return exited;
            });
        }

        bool find(const Graph &forbidden, SubgraphCallback callback) override {
            return graph.for_all_edges([&](VertexPair uv) {
                if (forbidden.has_edge(uv)) return false;
                Vertex u = uv.u;
                Vertex v = uv.v;

                w_candidates = graph.adj[u] & ~graph.adj[v] & ~forbidden.adj[u] & ~forbidden.adj[v];
                w_candidates[v] = false;

                bool exited = Graph::iterate(w_candidates, [&](Vertex w) {
                    return callback(Subgraph{u, v, w});
                });
                if (exited) return true;

                w_candidates = ~graph.adj[u] & graph.adj[v] & ~forbidden.adj[u] & ~forbidden.adj[v];
                w_candidates[u] = false;

                exited = Graph::iterate(w_candidates, [&](Vertex w) {
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
            Vertex u = uv.u;
            Vertex v = uv.v;

            if (graph.has_edge(uv)) {
                w_candidates = graph.adj[u] & ~graph.adj[v] & ~forbidden.adj[u] & ~forbidden.adj[v];
                w_candidates[v] = false;

                bool exited = Graph::iterate(w_candidates, [&](Vertex w) {
                    return callback(Subgraph{u, v, w});
                });
                if (exited) return true;

                w_candidates = ~graph.adj[u] & graph.adj[v] & ~forbidden.adj[u] & ~forbidden.adj[v];
                w_candidates[u] = false;

                exited = Graph::iterate(w_candidates, [&](Vertex w) {
                    return callback(Subgraph{v, u, w});
                });

                return exited;
            } else {
                w_candidates = graph.adj[u] & graph.adj[v] & ~forbidden.adj[u] & ~forbidden.adj[v];
                bool exited = Graph::iterate(w_candidates, [&](Vertex w) {
                    return callback(Subgraph{w, u, v});
                });

                return exited;
            }
        }

    };
}


#endif //CONCEPT_CENTERP3_H
