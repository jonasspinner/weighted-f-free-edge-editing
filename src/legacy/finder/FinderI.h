#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FINDERI_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FINDERI_H


#include "../../graph/Graph.h"
#include "../graph/Subgraph.h"
#include "../../Configuration.h"


class FinderI {

public:
    // TODO: Maybe convert to (const Subgraph &). This allows using the same Subgraph for finding and only copy it if
    //       needed.
    using SubgraphCallback = std::function<bool(Subgraph &&)>;

    // TODO: std::function is a major performance bottleneck. Use tlx::Delegate or replace all finder and subgraph code
    //       with templated classes.
    using VertexPairCallback = std::function<bool(VertexPair)>;

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
    virtual bool find(const Graph &graph, SubgraphCallback callback) = 0;

    /**
     * Find all forbidden subgraphs. Subgraphs sharing an edge with forbidden are ignored. When a subgraph is found,
     * callback(subgraph) is called.
     *
     * @param graph
     * @param forbidden
     * @param callback
     * @return
     */
    virtual bool find(const Graph &graph, const Graph &forbidden, SubgraphCallback callback) = 0;

    /**
     * Find all forbidden subgraphs having uv as an vertex pair. Subgraphs sharing an edge with forbidden are ignored.
     * When a subgraph is found, callback(subgraph) is called.
     *
     * @param uv
     * @param graph
     * @param callback
     * @return Whether the callback returned true once
     */
    virtual bool find_near(VertexPair uv, const Graph &graph, SubgraphCallback callback) = 0;

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
    virtual bool find_near(VertexPair uv, const Graph &graph, const Graph &forbidden, SubgraphCallback callback) = 0;


    virtual bool find_with_duplicates(const Graph &/*graph*/, const SubgraphCallback &/*callback*/) {
        throw std::runtime_error("FinderI::find_with_duplicates is not implemented");
    };

    virtual bool find_with_duplicates(const Graph &/*graph*/, const Graph &/*forbidden*/,
                                      const SubgraphCallback &/*callback*/) {
        throw std::runtime_error("FinderI::find_with_duplicates is not implemented");
    };

    virtual bool find_near_with_duplicates(VertexPair /*uv*/, const Graph &/*graph*/, const Graph &/*forbidden*/,
                                           const SubgraphCallback &/*callback*/) {
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
     * Contract with find_with_duplicates:
     *      All vertex pairs listed by this method are not forbidden. Vertex pairs that are not listed are permitted to
     *      be forbidden in the output of find_with_duplicates.
     *
     *      std::vector<Subgraph> subgraphs;
     *      Finder.find_with_duplicates(graph, forbidden, [](auto subgraph) {
     *          subgraphs.push_back(subgraph);
     *          return false;
     *      });
     *      for (auto subgraph : subgraphs) {
     *          finder.for_all_conversionless_edits(subgraph, [](auto uv) {
     *              assert(!forbidden.has_edge(uv));
     *          });
     *      }
     *
     *  TODO: Evaluate moving logic into Subgraph class.
     *  Note(jonas): Ideally the logic would be handled in subgraph.vertexPairsWithoutConversions(). This is not
     *      dependent on the type of forbidden subgraph, but also the set of forbidden subgraphs. The skipping of
     *      conversion for F = {C_l, P_l} is only possible because it is always the same pair of vertices, the degree
     *      one vertices in the graph, or a fixed edge in the cycle. This relies on the ordering of vertices in the
     *      Subgraph class, i.e. S[0] and S[l-1] are the fixed pair of vertices which, when editied, leads to a
     *      conversion.
     *
     * @return
     */
    virtual bool for_all_conversionless_edits(const Subgraph &/*subgraph*/,
            const VertexPairCallback &/*callback*/) const {
        throw std::runtime_error("FinderI::for_all_conversionless_edits is not implemented");
    }

    [[nodiscard]] Cost calculate_min_cost(const Subgraph &subgraph, const VertexPairMap<bool> &marked,
            const VertexPairMap<Cost> &costs) const {
        Cost min_cost = invalid_cost;
        for_all_conversionless_edits(subgraph, [&](auto uv) {
            if (!marked[uv]) {
                min_cost = std::min(min_cost, costs[uv]);
            }
            return false;
        });
        return min_cost;
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
        return [&](Vertex u) -> const Graph::AdjRow & { return graph.adj(u); };
    }

    static inline auto neighbors(const Graph &graph, const Graph &forbidden) {
        assert(graph.size() == forbidden.size());
        return [&](Vertex u) {
            auto result = graph.adj(u);
            result -= forbidden.adj(u);
            return result;
        };
    }

    static inline auto non_neighbors(const Graph &graph) {
        return [&](Vertex u) {
            auto result = graph.adj(u);
            result.flip();
            result[u] = false;
            return result;
        };
    }

    static inline auto non_neighbors(const Graph &graph, const Graph &forbidden) {
        assert(graph.size() == forbidden.size());
        return [&](Vertex u) {
            auto result = graph.adj(u);
            result.flip();
            result -= forbidden.adj(u);
            result[u] = false;
            return result;
        };
    }

    static inline auto valid_edge(const Graph &graph) {
        return [&](VertexPair xy) { return graph.has_edge(xy); };
    }

    static inline auto valid_edge(const Graph &graph, const Graph &forbidden) {
        assert(graph.size() == forbidden.size());
        return [&](VertexPair xy) { return graph.has_edge(xy) && !forbidden.has_edge(xy); };
    }

    static inline auto valid_non_edge(const Graph &graph) {
        return [&](VertexPair xy) { return !graph.has_edge(xy); };
    }

    static inline auto valid_non_edge(const Graph &graph, const Graph &forbidden) {
        assert(graph.size() == forbidden.size());
        return [&](VertexPair xy) { return !graph.has_edge(xy) && !forbidden.has_edge(xy); };
    }
};

#endif //WEIGHTED_F_FREE_EDGE_EDITING_FINDERI_H
