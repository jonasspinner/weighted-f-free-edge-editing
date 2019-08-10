//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_FINDERI_H
#define CONCEPT_FINDERI_H


#include "../graph/Graph.h"
#include "../graph/Subgraph.h"


class FinderI {

public:
    using SubgraphCallback = std::function<bool(Subgraph)>;

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


    static const Graph::AdjRow &neighbors(Vertex u, const Graph &graph) {
        return graph.adj[u];
    }

    static Graph::AdjRow neighbors(Vertex u, const Graph &graph, const Graph &forbidden) {
        return graph.adj[u] & ~forbidden.adj[u];
    }

    static Graph::AdjRow non_neighbors(Vertex u, const Graph &graph) {
        return ~graph.adj[u];
    }

    static Graph::AdjRow non_neighbors(Vertex u, const Graph &graph, const Graph &forbidden) {
        return ~(graph.adj[u] | forbidden.adj[u]);
    }

};

#endif //CONCEPT_FINDERI_H
