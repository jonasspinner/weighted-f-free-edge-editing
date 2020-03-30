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
    using VertexPairCallBack = std::function<bool(VertexPair)>;

public:
    virtual ~FinderI() = default;

    /**
     * Find all forbidden subgraphs. When a subgraph is found, callback(subgraph) is called.
     *
     * Contract: Every subgraph is listed exactly once.
     *      Example: P_3 1-2-3-1 list one of
     *          1-2-3-1, 2-3-1-2, 3-1-2-3,
     *          1-3-2-1, 2-1-3-2, 3-2-1-3
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


    virtual bool find_with_duplicates(const Graph& /*graph*/, const Graph &/*forbidden*/,
            SubgraphCallback /*callback*/) {
        throw std::runtime_error("FinderI::find_with_duplicates is not implemented");
    };


    /**
     * Iterate over all unmarked vertex pairs. A fixed vertex pair can be ignored if its edit would turn the given
     * forbidden subgraph into another.
     *
     * Example: F = (C_4, P_4)
     *      P_4: a-b-c-d   -> all except {a, d}. +{a, d} would result in a C_4
     *      C_4: a-b-c-d-a -> all except {a, d}. -{a, d} would result in a P_4
     *          Every other edge of the C_4 would also work, but it is assumed that the function will be called for all
     *          4 rotations (1-2-3-4-1, 2-3-4-1-2, 3-4-1-2-3, 4-1-2-3-4). As {a, d} is fixed, every edge will be
     *          excluded once.
     *
     * @return
     */
    virtual bool for_all_x_vertex_pairs(const Subgraph& /*subgraph*/, const VertexPairMap<bool> &/*marked*/,
                                        VertexPairCallBack /*callback*/) {
        throw std::runtime_error("FinderI::for_all_x_vertex_pairs is not implemented");
    }


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
