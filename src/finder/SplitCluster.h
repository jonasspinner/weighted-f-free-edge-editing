//
// Created by jonas on 02.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SPLITCLUSTER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SPLITCLUSTER_H


namespace Finder {
    class SplitCluster : public FinderI {
    private:
        Graph::AdjRow V;
        Graph::AdjRow W;
        Graph::AdjRow X;
        Graph::AdjRow Y;
    public:
        explicit SplitCluster(const Graph &graph_ref) : FinderI(graph_ref), V(graph.size()), W(graph.size()), X(graph.size()), Y(graph.size()) {}

        bool find(SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return  graph.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair uv) { return  graph.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv); };

            return find(callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find(const Graph& forbidden, SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return graph.m_adj[u] & ~forbidden.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u] & ~forbidden.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair uv) { return  graph.has_edge(uv) && !forbidden.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv) && !forbidden.has_edge(uv); };

            return find(callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find_near(VertexPair /*uv*/, SubgraphCallback /*callback*/) override {
            assert(false);
            return false;
        }

        bool find_near(VertexPair /*uv*/, const Graph& /*forbidden*/, SubgraphCallback /*callback*/) override  {
            assert(false);
            return false;
        }

        template <typename F, typename G, typename H, typename I>
        bool find(const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {

            // C4, C5, P5, Necktie, Bowtie

            // Necktie: {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {3, 5}}
            // Bowtie: {{1, 2}, {2, 3}, {1, 3}, {3, 4}, {4, 5}, {3, 5}}
            for (Vertex u : graph.vertices()) {
                V = neighbors(u);
                for (Vertex v : Graph::iterate(V)) {

                    // C4 / C5 / P5 / Bowtie / Necktie
                    W = neighbors(v) & non_neighbors(u);
                    for (Vertex w : Graph::iterate(W)) {

                        // C4
                        X = neighbors(u) & neighbors(w) & non_neighbors(v);
                        for (Vertex x : Graph::iterate(X)) {
                            assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_edge({u, x})); assert(valid_edge({v, w})); assert(valid_non_edge({v, x})); assert(valid_edge({w, x}));
                            if (callback({u, v, w, x})) return true;
                        }

                        // C5 / P5
                        X = neighbors(w) & non_neighbors(u) & non_neighbors(v);
                        Y = neighbors(u) & non_neighbors(v) & non_neighbors(w);
                        for (Vertex y : Graph::iterate(Y)) {
                            for (Vertex x : Graph::iterate(X)) {
                                if (valid_edge({x, y})) {
                                    // C5
                                    assert(valid_edge({y, u})); assert(valid_non_edge({y, v})); assert(valid_non_edge({y, w})); assert(valid_edge({y, x}));
                                    assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_non_edge({u, x})); assert(valid_edge({v, w})); assert(valid_non_edge({v, x})); assert(valid_edge({w, x}));
                                    if (callback({y, u, v, w, x})) return true;
                                } else if (valid_non_edge({x, y})) {
                                    // P5
                                    assert(valid_edge({y, u})); assert(valid_non_edge({y, v})); assert(valid_non_edge({y, w})); assert(valid_non_edge({y, x}));
                                    assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_non_edge({u, x})); assert(valid_edge({v, w})); assert(valid_non_edge({v, x})); assert(valid_edge({w, x}));
                                    if (callback({y, u, v, w, x})) return true;
                                }
                            }
                        }

                        // Necktie
                        X = neighbors(w) & neighbors(v) & non_neighbors(u);
                        for (Vertex x : Graph::iterate(X)) {
                            Y = neighbors(u) & non_neighbors(v) & non_neighbors(w) & non_neighbors(x);
                            for (Vertex y : Graph::iterate(Y)) {
                                assert(valid_edge({y, u})); assert(valid_non_edge({y, v})); assert(valid_non_edge({y, w})); assert(valid_non_edge({y, x}));
                                assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_non_edge({u, x})); assert(valid_edge({v, w})); assert(valid_edge({v, x})); assert(valid_edge({w, x}));
                                if (callback({y, u, v, w, x})) return true;
                            }
                        }

                        X = neighbors(w) & non_neighbors(v) & non_neighbors(u);
                        for (Vertex x : Graph::iterate(X)) {
                            Y = neighbors(u) & neighbors(v) & non_neighbors(w) & non_neighbors(x);
                            for (Vertex y : Graph::iterate(Y)) {
                                assert(valid_edge({y, u})); assert(valid_edge({y, v})); assert(valid_non_edge({y, w})); assert(valid_non_edge({y, x}));
                                assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_non_edge({u, x})); assert(valid_edge({v, w})); assert(valid_non_edge({v, x})); assert(valid_edge({w, x}));
                                if (callback({y, u, v, w, x})) return true;
                            }
                        }

                        // Bowtie
                        X = neighbors(w) & neighbors(v) & non_neighbors(u);
                        Y = neighbors(u) & neighbors(v) & non_neighbors(w);
                        for (Vertex x : Graph::iterate(X)) {
                            for (Vertex y : Graph::iterate(Y)) {
                                if (!valid_non_edge({x, y})) continue;
                                assert(valid_edge({y, u})); assert(valid_edge({y, v})); assert(valid_non_edge({y, w})); assert(valid_non_edge({y, x}));
                                assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_non_edge({u, x})); assert(valid_edge({v, w})); assert(valid_edge({v, x})); assert(valid_edge({w, x}));
                                if (callback({y, u, v, w, x})) return true;
                            }
                        }
                    }
                }
            }
            return false;
        }
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SPLITCLUSTER_H
