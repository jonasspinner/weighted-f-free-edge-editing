//
// Created by jonas on 03.07.19.
//

#ifndef CONCEPT_CENTERC4P4_H
#define CONCEPT_CENTERC4P4_H


#include "../graph/Graph.h"
#include "../graph/Subgraph.h"

namespace Finder {
    class CenterC4P4 : public FinderI {
        Graph::AdjRow V;
        Graph::AdjRow A;
        Graph::AdjRow B;

    public:
        explicit CenterC4P4(const Graph &graph_ref) : FinderI(graph_ref), V(graph.size()), A(graph.size()), B(graph.size()) {}

        template<typename F, typename G, typename H, typename I>
        bool find(const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {

            for (Vertex u : graph.vertices()) {
                V = neighbors(u);
                for (Vertex v : Graph::iterate(V)) {

                    A = neighbors(u) & non_neighbors(v);
                    B = neighbors(v) & non_neighbors(u);

                    for (Vertex a : Graph::iterate(A)) {
                        for (Vertex b : Graph::iterate(B)) {

                            assert(valid_edge({a, u}));
                            assert(valid_non_edge({a, v}));
                            // assert(valid({a, b}));
                            assert(valid_edge({u, v}));
                            assert(valid_non_edge({u, b}));
                            assert(valid_edge({v, b}));

                            if (valid_non_edge({a, b}) && a < b) {
                                // P_4
                                if (callback(Subgraph{a, u, v, b})) return true;
                            }
                            if (valid_edge({a, b}) && a < std::min({u, v, b}) && u < b) {
                                // C_4
                                if (callback(Subgraph{a, u, v, b})) return true;
                            }
                        }
                    }
                }
            }
            return false;
        }

        bool find(SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return  graph.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair uv) { return  graph.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv); };

            return find(callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find(const Graph &forbidden, SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return  graph.m_adj[u] - forbidden.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u] - forbidden.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair uv) { return  graph.has_edge(uv) && !forbidden.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv) && !forbidden.has_edge(uv); };

            return find(callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        template<typename F, typename G, typename H, typename I>
        bool find_near(VertexPair uv, const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge,
                       I valid_non_edge) {

            Vertex u = uv.u, v = uv.v;
            assert(u < v);

            if (valid_non_edge(uv)) {
                // case P_4 u, v are at positions 1 and 4 in the path
                A = neighbors(u) & non_neighbors(v);
                B = neighbors(v) & non_neighbors(u);

                for (Vertex a : Graph::iterate(A)) {
                    for (Vertex b : Graph::iterate(B)) {
                        if (valid_edge({a, b})) {
                            assert(valid_edge({u, a}));
                            assert(valid_non_edge({u, b}));
                            assert(valid_non_edge({u, v}));
                            assert(valid_edge({a, b}));
                            assert(valid_non_edge({a, v}));
                            assert(valid_edge({b, v}));

                            // P_4
                            if (callback(Subgraph{u, a, b, v})) return true;
                        }
                    }
                }

                // case P_4, C_4 u, v are at positions 1 and 3 in the path
                A = neighbors(u) & neighbors(v);
                for (Vertex a : Graph::iterate(A)) {
                    B = neighbors(v) & non_neighbors(a);
                    for (Vertex b : Graph::iterate(B)) {

                        assert(valid_edge({u, a}));
                        assert(valid_non_edge({u, v}));
                        // assert(valid({u, b}));
                        assert(valid_edge({a, v}));
                        assert(valid_non_edge({a, b}));
                        assert(valid_edge({v, b}));

                        if (valid_edge({u, b}) && u < std::min({a, v, b})) {
                            // C_4
                            if (callback(Subgraph{u, a, v, b})) return true;
                        } else if (valid_non_edge({u, b})) {
                            // P_4
                            if (callback(Subgraph{u, a, v, b})) return true;
                        }
                    }
                }

                // case P_4, C_4 u, v are at positions 2 and 4 in the path
                B = neighbors(u) & neighbors(v);
                for (Vertex b : Graph::iterate(B)) {
                    A = neighbors(u) & non_neighbors(b);
                    for (Vertex a : Graph::iterate(A)) {

                        assert(valid_edge({a, u}));
                        assert(valid_non_edge({a, b}));
                        // assert(valid({a, v}));
                        assert(valid_edge({u, b}));
                        assert(valid_non_edge({u, v}));
                        assert(valid_edge({b, v}));

                        if (valid_edge({a, v}) && a < std::min({u, b, v})) {
                            // C_4
                            if (callback(Subgraph{a, u, b, v})) return true;
                        } else if (valid_non_edge({a, v})) {
                            // P_4
                            if (callback(Subgraph{a, u, b, v})) return true;
                        }
                    }
                }

            } else if (valid_edge(uv)) {
                // case C_4 u, v are at positions 1 and 4 in the path
                A = neighbors(u) & non_neighbors(v);
                B = neighbors(v) & non_neighbors(u);

                for (Vertex a : Graph::iterate(A)) {
                    for (Vertex b : Graph::iterate(B)) {
                        if (valid_edge({a, b})) {

                            assert(valid_edge({u, a}));
                            assert(valid_non_edge({u, b}));
                            assert(valid_edge({u, v}));
                            assert(valid_edge({a, b}));
                            assert(valid_non_edge({a, v}));
                            assert(valid_edge({b, v}));

                            // C_4
                            if (callback(Subgraph{u, a, b, v})) return true;
                        }
                    }
                }

                // case C_4, P_4 u, v are at positions 2 and 3 in the path
                A = neighbors(u) & non_neighbors(v);
                B = neighbors(v) & non_neighbors(u);

                for (Vertex a : Graph::iterate(A)) {
                    for (Vertex b : Graph::iterate(B)) {

                        assert(valid_edge({a, u}));
                        assert(valid_non_edge({a, v}));
                        // assert(valid({a, b}));
                        assert(valid_edge({u, v}));
                        assert(valid_non_edge({u, b}));
                        assert(valid_edge({v, b}));

                        if (valid_edge({a, b}) && a < std::min({u, v, b})) {
                            // C_4
                            if (callback(Subgraph{a, u, v, b})) return true;
                        } else if (valid_non_edge({a, b})) {
                            // P_4
                            if (callback(Subgraph{a, u, v, b})) return true;
                        }
                    }
                }

                // case C_4, P_4 u, v are at positions 1 and 2 in the path
                A = neighbors(v) & non_neighbors(u);

                for (Vertex a : Graph::iterate(A)) {
                    B = neighbors(a) & non_neighbors(v);
                    B[u] = false;

                    for (Vertex b : Graph::iterate(B)) {
                        assert(valid_edge({u, v}));
                        assert(valid_non_edge({u, a}));
                        // assert(valid({u, b}));
                        assert(valid_edge({v, a}));
                        assert(valid_non_edge({v, b}));
                        assert(valid_edge({a, b}));

                        if (valid_edge({u, b}) && u < std::min({v, a, b})) {
                            // C_4
                            if (callback(Subgraph{u, v, a, b})) return true;
                        } else if (valid_non_edge({u, b})) {
                            // P_4
                            if (callback(Subgraph{u, v, a, b})) return true;
                        }
                    }
                }

                // case C_4, P_4 u, v are at positions 3 and 4 in the path
                A = neighbors(u) & non_neighbors(v);

                for (Vertex a : Graph::iterate(A)) {
                    B = neighbors(a) & non_neighbors(u);
                    B[v] = false;

                    for (Vertex b : Graph::iterate(B)) {
                        assert(valid_edge({b, a}));
                        assert(valid_non_edge({b, u}));
                        // assert(valid({b, v}));
                        assert(valid_edge({a, u}));
                        assert(valid_non_edge({a, v}));
                        assert(valid_edge({u, v}));

                        if (valid_edge({b, v}) && b < std::min({a, u, v})) {
                            // C_4
                            if (callback(Subgraph{b, a, u, v})) return true;
                        } else if (valid_non_edge({b, v})) {
                            // P_4
                            if (callback(Subgraph{b, a, u, v})) return true;
                        }
                    }
                }
            }
            return false;
        }

        bool find_near(VertexPair uv, SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return  graph.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair xy) { return  graph.has_edge(xy); };
            auto valid_non_edge = [&](VertexPair xy) { return !graph.has_edge(xy); };

            return find_near(uv, [&](Subgraph &&subgraph) {
                verify(subgraph, graph);
                return callback(std::move(subgraph));
            }, neighbors, non_neighbors, valid_edge, valid_non_edge);
        };

        bool find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return graph.m_adj[u] - forbidden.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u] - forbidden.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair xy) { return  graph.has_edge(xy) && !forbidden.has_edge(xy); };
            auto valid_non_edge = [&](VertexPair xy) { return !graph.has_edge(xy) && !forbidden.has_edge(xy); };

            return find_near(uv, [&](Subgraph &&subgraph) {
                verify(subgraph, graph);
                return callback(std::move(subgraph));
            }, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }
    };
}


#endif //CONCEPT_CENTERC4P4_H
