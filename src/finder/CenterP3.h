//
// Created by jonas on 23.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_CENTERP3_H
#define WEIGHTED_F_FREE_EDGE_EDITING_CENTERP3_H


#include "../interfaces/FinderI.h"


namespace Finder {
    class CenterP3 : public FinderI {

        Graph::AdjRow W;

    public:
        explicit CenterP3(const Graph &graph_ref) : FinderI(graph_ref), W(graph.size()) {}

        bool find(SubgraphCallback callback) override {
            auto neighbors =      [&](Vertex u)      { return  graph.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair xy) { return  graph.has_edge(xy); };
            auto valid_non_edge = [&](VertexPair xy) { return !graph.has_edge(xy); };

            return find(callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find(const Graph &forbidden, SubgraphCallback callback) override {
            auto neighbors =      [&](Vertex u)      { return  graph.m_adj[u] - forbidden.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u] - forbidden.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair xy) { return  graph.has_edge(xy) && !forbidden.has_edge(xy); };
            auto valid_non_edge = [&](VertexPair xy) { return !graph.has_edge(xy) && !forbidden.has_edge(xy); };

            return find(callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find_near(VertexPair uv, SubgraphCallback callback) override {
            auto neighbors =      [&](Vertex u)      { return  graph.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair xy) { return  graph.has_edge(xy); };
            auto valid_non_edge = [&](VertexPair xy) { return !graph.has_edge(xy); };

            return find_near(uv, callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) override {
            auto neighbors =      [&](Vertex u)      { return  graph.m_adj[u] - forbidden.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u] - forbidden.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair xy) { return  graph.has_edge(xy) && !forbidden.has_edge(xy); };
            auto valid_non_edge = [&](VertexPair xy) { return !graph.has_edge(xy) && !forbidden.has_edge(xy); };

            return find_near(uv, callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

    private:
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
        /**
         * Find paths of length 3. Each path is listed exactly once.
         */
        template<typename F, typename G, typename H, typename I>
        bool find(const SubgraphCallback& callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
#pragma GCC diagnostic pop
            for (Vertex u : graph.vertices()) {
                for (Vertex v : graph.neighbors(u)) {
                    if (!valid_edge({u, v})) continue;

                    W = neighbors(v) & non_neighbors(u);
                    for (Vertex w : Graph::iterate(W)) {
                        assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_edge({v, w}));
                        if (u < w) {
                            if (callback(Subgraph{u, v, w})) return true;
                        }
                    }
                }
            }
            return false;
        }

        template<typename F, typename G, typename H, typename I>
        bool find_near(VertexPair uv, const SubgraphCallback& callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {

            Vertex u = uv.u;
            Vertex v = uv.v;

            if (valid_edge(uv)) {
                W = neighbors(u) & non_neighbors(v);

                for (Vertex w : Graph::iterate(W)) {
                    assert(valid_edge({w, u})); assert(valid_non_edge({w, v})); assert(valid_edge({u, v}));
                    if (callback(Subgraph{w, u, v})) return true;
                }

                W = neighbors(v) & non_neighbors(u);

                for (Vertex w : Graph::iterate(W)) {
                    assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_edge({v, w}));
                    if (callback(Subgraph{u, v, w})) return true;
                }

            } else if (valid_non_edge(uv)) {

                W = neighbors(u) & neighbors(v);
                for (Vertex w : Graph::iterate(W)) {
                    assert(valid_edge({u, w})); assert(valid_non_edge({u, v})); assert(valid_edge({w, v}));
                    if (callback(Subgraph{u, w, v})) return true;
                }
            }
            return false;
        }

    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_CENTERP3_H
