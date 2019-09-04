//
// Created by jonas on 02.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINT_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINT_H


namespace Finder {
    template <size_t length>
    class Endpoint : public FinderI {
    public:
        explicit Endpoint(const Graph &graph_ref) : FinderI(graph_ref) {};

        bool find(SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return  graph.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair uv) { return  graph.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv); };

            return detail::FindImpl<length, (length > 3)>::find(graph, callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find(const Graph& forbidden, SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return graph.m_adj[u] & ~forbidden.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u] & ~forbidden.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair uv) { return  graph.has_edge(uv) && !forbidden.has_edge(uv); };
            auto valid_non_edge = [&](VertexPair uv) { return !graph.has_edge(uv) && !forbidden.has_edge(uv); };

            return detail::FindImpl<length, (length > 3)>::find(graph, callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find_near(VertexPair uv, SubgraphCallback callback) override {

            auto neighbors =      [&](Vertex u)      { return graph.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair xy) { return  graph.has_edge(xy); };
            auto valid_non_edge = [&](VertexPair xy) { return !graph.has_edge(xy); };

            return detail::FindNearImpl<length, (length > 3)>::find_near(graph, uv, callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        bool find_near(VertexPair uv, const Graph& forbidden, SubgraphCallback callback) override  {

            auto neighbors =      [&](Vertex u)      { return graph.m_adj[u] & ~forbidden.m_adj[u]; };
            auto non_neighbors =  [&](Vertex u)      { auto result = ~graph.m_adj[u] & ~forbidden.m_adj[u]; result[u] = false; return result; };
            auto valid_edge =     [&](VertexPair xy) { return  graph.has_edge(xy) && !forbidden.has_edge(xy); };
            auto valid_non_edge = [&](VertexPair xy) { return !graph.has_edge(xy) && !forbidden.has_edge(xy); };

            return detail::FindNearImpl<length, (length > 3)>::find_near(graph, uv, callback, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINT_H
