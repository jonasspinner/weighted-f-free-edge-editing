//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_FINDERI_H
#define CONCEPT_FINDERI_H


#include "../graph/Graph.h"
#include "../graph/Subgraph.h"


class FinderI {

public:
    using SubgraphCallback = std::function<bool(Subgraph&&)>;

protected:
    const Graph &graph;

public:
    explicit FinderI(const Graph &graph_ref) : graph(graph_ref) {}

    virtual ~FinderI() = default;

    /**
     * Find all forbidden subgraphs. When a subgraph is found, callback(subgraph) is called.
     *
     * @param callback
     * @return Whether the callback returned true once
     */
    virtual bool find(SubgraphCallback callback) = 0;

    /**
     * Find all forbidden subgraphs. Subgraphs sharing an edge with forbidden are ignored. When a subgraph is found,
     * callback(subgraph) is called.
     *
     * @param forbidden
     * @param callback
     * @return
     */
    virtual bool find(const Graph &forbidden, SubgraphCallback callback) = 0;

    /**
     * Find all forbidden subgraphs having uv as an vertex pair. Subgraphs sharing an edge with forbidden are ignored.
     * When a subgraph is found, callback(subgraph) is called.
     *
     * @param uv
     * @param callback
     * @return Whether the callback returned true once
     */
    virtual bool find_near(VertexPair uv, SubgraphCallback callback) = 0;

    /**
     * Find all forbidden subgraphs having uv as an vertex pair. Subgraphs sharing an edge with forbidden are ignored.
     * When a subgraph is found, callback(subgraph) is called.
     *
     * @param uv
     * @param forbidden
     * @param callback
     * @return Whether the callback returned true once
     */
    virtual bool find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) = 0;

protected:
    /*
     * Utility functions for generating lambda functions for finder templates.
     */

    static inline auto neighbors(const Graph &graph) {
        return [&](Vertex u) { return  graph.m_adj[u]; };
    }

    static inline auto neighbors(const Graph &graph, const Graph &forbidden) {
        return [&](Vertex u) { return graph.m_adj[u] - forbidden.m_adj[u]; };
    }

    static inline auto non_neighbors(const Graph &graph) {
        return [&](Vertex u) { auto result = ~graph.m_adj[u]; result[u] = false; return result; };
    }

    static inline auto non_neighbors(const Graph &graph, const Graph &forbidden) {
        return [&](Vertex u) { auto result = ~graph.m_adj[u] - forbidden.m_adj[u]; result[u] = false; return result; };
    }

    static inline auto valid_edge(const Graph &graph) {
        return [&](VertexPair xy) { return graph.has_edge(xy); };
    }

    static inline auto valid_edge(const Graph &graph, const Graph &forbidden) {
        return [&](VertexPair xy) { return graph.has_edge(xy) && !forbidden.has_edge(xy); };
    }

    static inline auto valid_non_edge(const Graph &graph) {
        return [&](VertexPair xy) { return !graph.has_edge(xy); };
    }

    static inline auto valid_non_edge(const Graph &graph, const Graph &forbidden) {
        return [&](VertexPair xy) { return !graph.has_edge(xy) && !forbidden.has_edge(xy); };
    }
};

#endif //CONCEPT_FINDERI_H
