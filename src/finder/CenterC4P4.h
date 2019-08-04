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
        explicit CenterC4P4(const Graph& graph) : FinderI(graph), w_candidate(graph.size()), x_candidate(
                graph.size()) {}

        bool find(SubgraphCallback callback) override {
            /* w?x
             * | |
             * u-v
             */
            /*
            return graph.for_all_edges([&](VertexPair uv) {
                Vertex u = uv.u;
                Vertex v = uv.v;

                // Adj(u) \ Adj(v) - {v}
                // auto w_candidates = graph.adj[u] & ~graph.adj[v];

                // Adj(v) \ Adj(u) - {u}
                // auto x_candidates = graph.adj[v] & ~graph.adj[u];

                w_candidate = graph.adj[u] & ~graph.adj[v];
                w_candidate[v] = false;

                x_candidate = graph.adj[v] & ~graph.adj[u];
                x_candidate[u] = false;

                return Graph::iterate(w_candidate, [&](Vertex w) {
                    return Graph::iterate(x_candidate, [&](Vertex x) {
                        if (w == x) return false;
                        if (u > v || w > x) return false;
                        // if (graph.has_edge({w, x}) && (u > v || u > w || u > x)) return false;
                        return callback(Subgraph{u, v, w, x});
                    });
                });
            });*/
            return graph.for_all_vertices([&](Vertex p_1) {
                return graph.for_neighbors_of(p_1, [&](Vertex p_2) {
                    if (p_1 >= p_2) return false; // TODO: not specified in paper
                    auto A = graph.adj[p_1] & ~graph.adj[p_2];
                    //A[p_2] = false;
                    auto B = graph.adj[p_2] & ~graph.adj[p_1];
                    //B[p_1] = false;
                    return Graph::iterate(A, [&](Vertex a) {
                        return Graph::iterate(B, [&](Vertex b){
                            if (!graph.has_edge({a, b})) return callback(Subgraph{a, p_1, p_2, b});
                            if (graph.has_edge({a, b}) && (p_1 < p_2) && (p_1 < a) && (p_1 < b) && (p_2 < a)) return callback(Subgraph{p_1, p_2, b, a});
                            return false;
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


#endif //CONCEPT_CENTERC4P4_H
