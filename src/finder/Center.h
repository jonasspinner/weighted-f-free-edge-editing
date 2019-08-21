//
// Created by jonas on 31.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H


#include "../interfaces/FinderI.h"
#include "center/find.h"
#include "center/find_near.h"


namespace detail {

    template <int length>
    class Center : public FinderI {
        static_assert(length > 1);

    public:
        explicit Center(const Graph &graph_ref) : FinderI(graph_ref) {}

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

using CenterRecC4P4 = detail::Center<4>;
using CenterRecP3 = detail::Center<3>;

#endif //WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H
