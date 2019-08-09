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
            /*
            return graph.for_all_vertices([&](Vertex p_1) {
                return graph.for_neighbors_of(p_1, [&](Vertex p_2) {

                    auto A = neighbors(p_1) & non_neighbors(p_2); A[p_2] = false;
                    auto B = neighbors(p_2) & non_neighbors(p_1); B[p_1] = false;

                    return Graph::iterate(A, B, [&](Vertex a, Vertex b) {
                        if (valid_non_edge({a, b}) && a < b) {
                            assert(valid_edge({a, p_1})); assert(valid_non_edge({a, p_2})); assert(valid_non_edge({a, b})); assert(valid_edge({p_1, p_2})); assert(valid_non_edge({p_1, b})); assert(valid_edge({p_2, b}));
                            return callback(Subgraph{a, p_1, p_2, b});
                        }
                        if (valid_edge({a, b}) && (p_1 < p_2) && (p_1 < b) && (p_1 < a) && (p_2 < a)) {
                            assert(valid_edge({a, p_1})); assert(valid_non_edge({a, p_2})); assert(valid_edge({a, b})); assert(valid_edge({p_1, p_2})); assert(valid_non_edge({p_1, b})); assert(valid_edge({p_2, b}));
                            return callback(Subgraph{p_1, p_2, b, a});
                        }
                        return false;
                    });
                });
            });*/

            for (Vertex u : graph.vertices()) {
                for (Vertex v : graph.neighbors(u)) {

                    auto A = neighbors(u) & non_neighbors(v);
                    auto B = neighbors(v) & non_neighbors(u);

                    for (Vertex a : Graph::vertices(A)) {
                        for (Vertex b : Graph::vertices(B)) {

                            if (valid_non_edge({a, b}) && a < b) {
                                assert(valid_edge({a, u})); assert(valid_non_edge({a, v})); assert(valid_non_edge({a, b})); assert(valid_edge({u, v})); assert(valid_non_edge({u, b})); assert(valid_edge({v, b}));
                                if(callback(Subgraph{a, u, v, b})) return true;
                            }
                            if (valid_edge({a, b}) && (u < v) && (u < b) && (u < a) && (v < a)) {
                                assert(valid_edge({u, v})); assert(valid_non_edge({u, b})); assert(valid_edge({u, a})); assert(valid_edge({v, b})); assert(valid_non_edge({v, a})); assert(valid_edge({b, a}));
                                if(callback(Subgraph{u, v, b, a})) return true;
                            }
                        }
                    }
                }
            }
            return false;
        }

        bool find(SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return  graph.adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair uv) { return  graph.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv); };

            return find(callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find(const Graph &forbidden, SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return  graph.adj[u] & ~forbidden.adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.adj[u] & ~forbidden.adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair uv) { return  graph.has_edge(uv) && !forbidden.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv) && !forbidden.has_edge(uv); };

            return find(callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        template <typename F, typename G, typename H, typename I>
        bool find_near(VertexPair uv, const SubgraphCallback& callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            Vertex u = uv.u, v = uv.v;
            if (valid_non_edge(uv)) {
                // case P_4
                auto A = neighbors(u) & non_neighbors(v);
                auto B = neighbors(v) & non_neighbors(u);
                return Graph::iterate(A, B, [&](Vertex a, Vertex b) {
                    if (valid_edge({a, b})) {
                        assert(valid_edge({u, a})); assert(valid_non_edge({u, b})); assert(valid_non_edge({u, v})); assert(valid_edge({a, b})); assert(valid_non_edge({a, v})); assert(valid_edge({b, v}));
                        return callback(Subgraph{u, a, b, v});
                    }
                    return false;
                });
            } else if (valid_edge(uv)){
                // case C_4, P_4 uv is in the middle
                auto A = neighbors(u) & non_neighbors(v);
                auto B = neighbors(v) & non_neighbors(u);
                bool exit_early = Graph::iterate(A, B, [&](Vertex a, Vertex b) {
                    assert(valid_edge({a, u})); assert(valid_non_edge({a, v})); assert(valid_edge({u, v})); assert(valid_non_edge({u, b})); assert(valid_edge({v, b}));
                    return callback(Subgraph{a, u, v, b});
                });
                if (exit_early) return true;

                // case C_4, P_4 uv is on the left
                w_candidate = neighbors(v) & non_neighbors(u);
                exit_early = Graph::iterate(w_candidate, [&](Vertex w) {
                    x_candidate = neighbors(w) & non_neighbors(v); x_candidate[u] = false;
                    return Graph::iterate(x_candidate, [&](Vertex x) {
                        assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_edge({v, w})); assert(valid_non_edge({v, x})); assert(valid_edge({w, x}));
                        return callback(Subgraph{u, v, w, x});
                    });
                });
                if (exit_early) return true;

                // case C_4, P_4 uv is on the right
                x_candidate = neighbors(u) & non_neighbors(v);
                exit_early = Graph::iterate(x_candidate, [&](Vertex x) {
                    w_candidate = neighbors(x) & non_neighbors(u); w_candidate[v] = false;
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

            auto neighbors =      [&](Vertex u)      { return  graph.adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair uv) { return  graph.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv); };

            return find_near(uv, [&](Subgraph && subgraph) { verify(subgraph, graph); return callback(std::move(subgraph)); }, neighbors, non_neighbors, valid_edge, valid_non_edge);
        };

        bool find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return  graph.adj[u] & ~forbidden.adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.adj[u] & ~forbidden.adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair uv) { return  graph.has_edge(uv) && !forbidden.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv) && !forbidden.has_edge(uv); };

            return find_near(uv, [&](Subgraph && subgraph) { verify(subgraph, graph); return callback(std::move(subgraph)); }, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }
    };
}


#endif //CONCEPT_CENTERC4P4_H
