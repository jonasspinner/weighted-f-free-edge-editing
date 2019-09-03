//
// Created by jonas on 02.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SPLITGRAPH_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SPLITGRAPH_H


namespace Finder {
    class SplitGraph : public FinderI {
    private:
        Graph::AdjRow V;
        Graph::AdjRow W;
        Graph::AdjRow X;
        Graph::AdjRow Y;
    public:
        explicit SplitGraph(const Graph &graph_ref) : FinderI(graph_ref), V(graph.size()), W(graph.size()), X(graph.size()), Y(graph.size()) {}

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

            for (Vertex u : graph.vertices()) {
                V = neighbors(u);
                for (Vertex v : Graph::iterate(V)) {

                    // 2K2
                    W = non_neighbors(u) & non_neighbors(v);
                    for (Vertex w : Graph::iterate(W)) {
                        X = neighbors(w) & non_neighbors(u) & non_neighbors(v);
                        for (Vertex x : Graph::iterate(X)) {
                            assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_non_edge({u, x})); assert(valid_non_edge({v, w})); assert(valid_non_edge({v, x})); assert(valid_edge({w, x}));
                            if (callback({u, v, w, x})) return true;
                        }
                    }

                    // C4 / C5
                    W = neighbors(v) & non_neighbors(u);
                    Y = neighbors(u) & non_neighbors(v);
                    for (Vertex w : Graph::iterate(W)) {

                        // C4
                        X = neighbors(u) & neighbors(w) & non_neighbors(v);
                        for (Vertex x : Graph::iterate(X)) {
                            assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_edge({u, x})); assert(valid_edge({v, w})); assert(valid_non_edge({v, x})); assert(valid_edge({w, x}));
                            if (callback({u, v, w, x})) return true;
                        }

                        // C5
                        for (Vertex y : Graph::iterate(Y)) {
                            if (!valid_non_edge({y, w})) continue;
                            X = neighbors(y) & neighbors(w) & non_neighbors(u) & non_neighbors(v);
                            for (Vertex x : Graph::iterate(X)) {
                                assert(valid_edge({y, u})); assert(valid_non_edge({y, v})); assert(valid_non_edge({y, w})); assert(valid_edge({y, x}));
                                assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_non_edge({u, x})); assert(valid_edge({v, w})); assert(valid_non_edge({v, x})); assert(valid_edge({w, x}));
                                if (callback({y, u, v, w, x})) return true;
                            }
                        }
                    }
                }
            }
            return false;
        }

        template <typename F, typename G, typename H, typename I>
        bool find_near(VertexPair uv, const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            auto [u, v] = uv;

            if (valid_edge(uv)) {
                // 2K2
                W = non_neighbors(u) & non_neighbors(v);
                for (Vertex w : Graph::iterate(W)) {
                    X = neighbors(w) & non_neighbors(u) & non_neighbors(v);
                    for (Vertex x : Graph::iterate(X)) {
                        assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_non_edge({u, x})); assert(valid_non_edge({v, w})); assert(valid_non_edge({v, x})); assert(valid_edge({w, x}));
                        if (callback({u, v, w, x})) return true;
                    }
                }

                assert(false);


            } else if (valid_non_edge(uv)) {

            }
            return false;
        }
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SPLITGRAPH_H