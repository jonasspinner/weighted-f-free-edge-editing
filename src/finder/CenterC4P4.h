//
// Created by jonas on 03.07.19.
//

#ifndef CONCEPT_CENTERC4P4_H
#define CONCEPT_CENTERC4P4_H


#include "../graph/Graph.h"
#include "../graph/Subgraph.h"

namespace Finder {
    class CenterC4P4 : public FinderI {
        Graph::AdjRow w_candidate;
        Graph::AdjRow x_candidate;

    public:
        explicit CenterC4P4(const Graph& graph) : FinderI(graph), w_candidate(graph.n_vertices()), x_candidate(graph.n_vertices()) {}

        bool find(SubgraphCallback callback) override {
            /* w?x
             * | |
             * u-v
             */
            return graph.for_all_edges([&](VertexPair uv) {
                Vertex u = uv.u;
                Vertex v = uv.v;

                // Adj(u) \ Adj(v) - {v}
                // auto w_candidates = graph.adj[u] & ~graph.adj[v];

                // Adj(v) \ Adj(u) - {u}
                // auto x_candidates = graph.adj[v] & ~graph.adj[u];

                w_candidate = graph.adj[u] & ~graph.adj[v];
                w_candidate[v] = false;

                return Graph::iterate(w_candidate, [&](Vertex w) {
                    x_candidate = graph.adj[v] & ~graph.adj[u];
                    x_candidate[u] = false;

                    return Graph::iterate(x_candidate, [&](Vertex x) {
                        // if (graph.has_edge({w, x}) && (u > v || u > w || u > x)) return false;
                        return callback(Subgraph{u, v, w, x});
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


#endif //CONCEPT_CENTERC4P4_H
