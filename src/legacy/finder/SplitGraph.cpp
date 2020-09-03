//
// Created by jonas on 02.09.19.
//


#include "SplitGraph.h"
#include "../Configuration.h"


namespace Finder {
    /**
     * Calls callback on 2K_2, C_4 and C_5 subgraphs.
     *
     * @param graph
     * @param callback
     * @return
     */
    bool SplitGraph::find(const Graph& graph, SubgraphCallback callback) {
        return find(graph, callback, neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
    }

    /**
     * Calls callback on 2K_2, C_4 and C_5 subgraphs. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
     *
     * @param graph
     * @param forbidden
     * @param callback
     * @return
     */
    bool SplitGraph::find(const Graph& graph, const Graph& forbidden, SubgraphCallback callback) {
        return find(graph, callback, neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    /**
     * Not implemented.
     *
     * @return
     */
    bool SplitGraph::find_near(VertexPair /*uv*/, const Graph& /*graph*/, SubgraphCallback /*callback*/) {
        throw std::runtime_error("SplitGraph does not support find_near");
        return false;
    }

    /**
     * Not implemented.
     *
     * @return
     */
    bool SplitGraph::find_near(VertexPair /*uv*/, const Graph& /*graph*/, const Graph& /*forbidden*/, SubgraphCallback /*callback*/)  {
        throw std::runtime_error("SplitGraph does not support find_near");
        return false;
    }

    void SplitGraph::to_yaml(YAML::Emitter &out) const {
        using namespace YAML;
        out << BeginMap;
        out << Key << "name" << Value << "SplitGraph";
        out << Key << "forbidden_subgraphs" << Value << Options::FSG::C4_C5_2K2;
        out << EndMap;
    }

    template <typename F, typename G, typename H, typename I>
    bool SplitGraph::find(const Graph& graph, const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
        // TODO: avoid listing a subgraph more than once

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
    bool SplitGraph::find_near(VertexPair uv, const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
        assert(false);
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
}
