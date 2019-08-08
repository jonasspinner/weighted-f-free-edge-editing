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

        template <typename F, typename G, typename H, typename I>
        bool find(const SubgraphCallback& callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            /* w?x
             * | |
             * u-v
             */
            /*
            return graph.for_all_edges([&](VertexPair uv) {
                Vertex u = uv.u;
                Vertex v = uv.v;

                // Adj(u) \ Adj(v) - {v}
                w_candidate = neighbors(u) & non_neighbors(v); w_candidate[v] = false;

                // Adj(v) \ Adj(u) - {u}
                x_candidate = neighbors(v) & non_neighbors(u); x_candidate[u] = false;

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

                    auto A = neighbors(p_1) & non_neighbors(p_2); A[p_2] = false;
                    auto B = neighbors(p_2) & non_neighbors(p_1); B[p_1] = false;

                    return Graph::iterate(A, B, [&](Vertex a, Vertex b) {
                        if (valid_non_edge({a, b})) {
                            assert(valid_edge({a, p_1})); assert(valid_non_edge({a, p_2})); assert(valid_non_edge({a, b})); assert(valid_edge({p_1, p_2})); assert(valid_non_edge({p_1, b})); assert(valid_edge({p_2, b}));
                            return callback(Subgraph{a, p_1, p_2, b});
                        }
                        if (valid_edge({a, b}) && (p_1 < p_2) && (p_1 < a) && (p_1 < b) && (p_2 < a)) {
                            assert(valid_edge({a, p_1})); assert(valid_non_edge({a, p_2})); assert(valid_edge({a, b})); assert(valid_edge({p_1, p_2})); assert(valid_non_edge({p_1, b})); assert(valid_edge({p_2, b}));
                            return callback(Subgraph{p_1, p_2, b, a});
                        }
                        return false;
                    });
                });
            });
        }

        bool find(SubgraphCallback callback) override {
            auto neighbors = [&](Vertex u) { return graph.adj[u]; };
            auto non_neighbors = [&](Vertex u) { return ~graph.adj[u]; };
            auto valid_edge = [&](VertexPair uv) { return graph.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return ~graph.has_edge(uv); };
            return find(callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find(const Graph &forbidden, SubgraphCallback callback) override {
            auto neighbors = [&](Vertex u) { return graph.adj[u] & ~forbidden.adj[u]; };
            auto non_neighbors = [&](Vertex u) { return ~graph.adj[u] & ~forbidden.adj[u]; };
            auto valid_edge = [&](VertexPair uv) { return graph.has_edge(uv) && !forbidden.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv) && !forbidden.has_edge(uv); };
            return find(callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        template <typename F, typename G, typename H, typename I>
        bool find_near(VertexPair uv, const SubgraphCallback& callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            Vertex u = uv.u, v = uv.v;
            if (valid_non_edge(uv)) {
                // case P_4
                auto A = neighbors(u) & non_neighbors(v); A[v] = false;
                auto B = neighbors(v) & non_neighbors(u); B[u] = false;
                return Graph::iterate(A, B, [&](Vertex a, Vertex b) {
                    if (valid_edge({a, b})) {
                        assert(valid_edge({u, a})); assert(valid_non_edge({u, b})); assert(valid_non_edge({u, v})); assert(valid_edge({a, b})); assert(valid_non_edge({a, v})); assert(valid_edge({b, v}));
                        return callback(Subgraph{u, a, b, v});
                    }
                    return false;
                });
            } else if (valid_edge(uv)){
                // case C_4, P_4 uv is in the middle
                auto A = neighbors(u) & non_neighbors(v); A[v] = false;
                auto B = neighbors(v) & non_neighbors(u); B[u] = false;
                bool exit_early = Graph::iterate(A, B, [&](Vertex a, Vertex b) {
                    assert(valid_edge({a, u})); assert(valid_non_edge({a, v})); assert(valid_edge({u, v})); assert(valid_non_edge({u, b})); assert(valid_edge({v, b}));
                    return callback(Subgraph{a, u, v, b});
                });
                if (exit_early) return true;

                // case C_4, P_4 uv is on the left
                w_candidate = neighbors(v) & non_neighbors(u); w_candidate[u] = false;
                exit_early = Graph::iterate(w_candidate, [&](Vertex w) {
                    x_candidate = neighbors(w) & non_neighbors(v); x_candidate[v] = false; x_candidate[u] = false;
                    return Graph::iterate(x_candidate, [&](Vertex x) {
                        assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_edge({v, w})); assert(valid_non_edge({v, x})); assert(valid_edge({w, x}));
                        return callback(Subgraph{u, v, w, x});
                    });
                });
                if (exit_early) return true;

                // case C_4, P_4 uv is on the right
                x_candidate = neighbors(u) & non_neighbors(v); x_candidate[v] = false;
                exit_early = Graph::iterate(x_candidate, [&](Vertex x) {
                    w_candidate = neighbors(x) & non_neighbors(u); w_candidate[u] = false; w_candidate[v] = false;
                    return Graph::iterate(w_candidate, [&](Vertex w) {
                        assert(valid_edge({w, x})); assert(valid_non_edge({w, u})); assert(valid_edge({x, u})); assert(valid_non_edge({x, v})); assert(valid_edge({u, v}));
                        return callback(Subgraph{w, x, u, v});
                    });
                });
                return exit_early;
            }
            return false;
        }

        bool find_near(VertexPair uv, SubgraphCallback callback) override {
            auto neighbors = [&](Vertex u) { return graph.adj[u]; };
            auto non_neighbors = [&](Vertex u) { return ~graph.adj[u]; };
            auto valid_edge = [&](VertexPair uv) { return graph.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv); };
            return find_near(uv, [&](Subgraph && subgraph) { verify(subgraph, graph); return callback(std::move(subgraph)); }, neighbors, non_neighbors, valid_edge, valid_non_edge);
        };

        bool find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) override {
            auto neighbors = [&](Vertex u) { return graph.adj[u] & ~forbidden.adj[u]; };
            auto non_neighbors = [&](Vertex u) { return ~graph.adj[u] & ~forbidden.adj[u]; };
            auto valid_edge = [&](VertexPair uv) { return graph.has_edge(uv) && !forbidden.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv) && !forbidden.has_edge(uv); };
            return find_near(uv, [&](Subgraph && subgraph) { verify(subgraph, graph); return callback(std::move(subgraph)); }, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }
    };
}


#endif //CONCEPT_CENTERC4P4_H
