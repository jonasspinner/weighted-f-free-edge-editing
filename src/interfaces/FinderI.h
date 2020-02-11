//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_FINDERI_H
#define CONCEPT_FINDERI_H


#include "../graph/Graph.h"
#include "../graph/Subgraph.h"
#include "../Configuration.h"


class FinderI {

public:
    using SubgraphCallback = std::function<bool(Subgraph&&)>;

public:
    virtual ~FinderI() = default;

    /**
     * Find all forbidden subgraphs. When a subgraph is found, callback(subgraph) is called.
     *
     * @param graph
     * @param callback
     * @return Whether the callback returned true once
     */
    virtual bool find(const Graph& graph, SubgraphCallback callback) = 0;

    /**
     * Find all forbidden subgraphs. Subgraphs sharing an edge with forbidden are ignored. When a subgraph is found,
     * callback(subgraph) is called.
     *
     * @param graph
     * @param forbidden
     * @param callback
     * @return
     */
    virtual bool find(const Graph& graph, const Graph &forbidden, SubgraphCallback callback) = 0;

    /**
     * Find all forbidden subgraphs having uv as an vertex pair. Subgraphs sharing an edge with forbidden are ignored.
     * When a subgraph is found, callback(subgraph) is called.
     *
     * @param uv
     * @param graph
     * @param callback
     * @return Whether the callback returned true once
     */
    virtual bool find_near(VertexPair uv, const Graph& graph, SubgraphCallback callback) = 0;

    /**
     * Find all forbidden subgraphs having uv as an vertex pair. Subgraphs sharing an edge with forbidden are ignored.
     * When a subgraph is found, callback(subgraph) is called.
     *
     * @param uv
     * @param graph
     * @param forbidden
     * @param callback
     * @return Whether the callback returned true once
     */
    virtual bool find_near(VertexPair uv, const Graph& graph, const Graph &forbidden, SubgraphCallback callback) = 0;

    [[nodiscard]] virtual Options::FSG forbidden_subgraphs() const = 0;

    [[nodiscard]] virtual std::string name() const = 0;

    virtual void to_yaml(YAML::Emitter &out) const = 0;

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const FinderI &finder) {
        finder.to_yaml(out);
        return out;
    }

protected:
    /*
     * Utility functions for generating lambda functions for finder templates.
     */

    static inline auto neighbors(const Graph &graph) {
        return [&](Vertex u) -> const Graph::AdjRow& { return  graph.m_adj[u]; };
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
        return [&](VertexPair xy) { return graph.hasEdge(xy); };
    }

    static inline auto valid_edge(const Graph &graph, const Graph &forbidden) {
        return [&](VertexPair xy) { return graph.hasEdge(xy) && !forbidden.hasEdge(xy); };
    }

    static inline auto valid_non_edge(const Graph &graph) {
        return [&](VertexPair xy) { return !graph.hasEdge(xy); };
    }

    static inline auto valid_non_edge(const Graph &graph, const Graph &forbidden) {
        return [&](VertexPair xy) { return !graph.hasEdge(xy) && !forbidden.hasEdge(xy); };
    }
};

#endif //CONCEPT_FINDERI_H
