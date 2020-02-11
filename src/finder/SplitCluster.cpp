//
// Created by jonas on 02.09.19.
//


#include "SplitCluster.h"
#include "../Configuration.h"


namespace Finder {

    /**
     * Calls callback on C4, C5, P5, Necktie and Bowtie subgraphs.
     *
     * @param graph
     * @param callback
     * @return
     */
    bool SplitCluster::find(const Graph& graph, SubgraphCallback callback) {
        return find(graph, callback, neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
    }

    /**
     * Calls callback on C4, C5, P5, Necktie and Bowtie subgraphs. Subgraphs sharing a vertex pair with the graph
     * forbidden are ignored.
     *
     * @param graph
     * @param forbidden
     * @param callback
     * @return
     */
    bool SplitCluster::find(const Graph& graph, const Graph& forbidden, SubgraphCallback callback) {
        return find(graph, callback, neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    /**
     * Not implemented.
     *
     * @return
     */
    bool SplitCluster::find_near(VertexPair /*uv*/, const Graph& /*graph*/, SubgraphCallback /*callback*/) {
        throw std::runtime_error("SplitCluster does not support find_near");
        return false;
    }

    /**
     * Not implemented.
     *
     * @return
     */
    bool SplitCluster::find_near(VertexPair /*uv*/, const Graph& /*graph*/, const Graph& /*forbidden*/, SubgraphCallback /*callback*/) {
        throw std::runtime_error("SplitCluster does not support find_near");
        return false;
    }

    void SplitCluster::to_yaml(YAML::Emitter &out) const {
        using namespace YAML;
        out << BeginMap;
        out << Key << "name" << Value << "SplitCluster";
        out << Key << "forbidden_subgraphs" << Value << Options::FSG::C4_C5_P5_Bowtie_Necktie;
        out << EndMap;
    }

    template <typename F, typename G, typename H, typename I>
    bool SplitCluster::find(const Graph& graph, const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {

        // C4, C5, P5, Necktie, Bowtie

        //  Necktie       Bowtie
        //   1---2         1---2
        //    \   \         \___\
        //     \---3            >3
        //        /         /¯¯ /
        //   5---4         5---4

        for (Vertex u : graph.vertices()) {
            V = neighbors(u);
            for (Vertex v : Graph::iterate(V)) {

                // C4 / C5 / P5 / Bowtie / Necktie
                W = neighbors(v) & non_neighbors(u);
                for (Vertex w : Graph::iterate(W)) {

                    // C4
                    X = neighbors(u) & neighbors(w) & non_neighbors(v);
                    for (Vertex x : Graph::iterate(X)) {
                        assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_edge({u, x})); assert(valid_edge({v, w})); assert(valid_non_edge({v, x})); assert(valid_edge({w, x}));

                        if (callback({u, v, w, x})) return true;
                    }

                    // C5 / P5
                    X = neighbors(w) & non_neighbors(u) & non_neighbors(v);
                    Y = neighbors(u) & non_neighbors(v) & non_neighbors(w);
                    for (Vertex y : Graph::iterate(Y)) {
                        for (Vertex x : Graph::iterate(X)) {
                            if (valid_edge({x, y})) {
                                // C5
                                assert(valid_edge({y, u})); assert(valid_non_edge({y, v})); assert(valid_non_edge({y, w})); assert(valid_edge({y, x}));
                                assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_non_edge({u, x})); assert(valid_edge({v, w})); assert(valid_non_edge({v, x})); assert(valid_edge({w, x}));
                                if (callback({y, u, v, w, x})) return true;
                            } else if (valid_non_edge({x, y})) {
                                // P5
                                assert(valid_edge({y, u})); assert(valid_non_edge({y, v})); assert(valid_non_edge({y, w})); assert(valid_non_edge({y, x}));
                                assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_non_edge({u, x})); assert(valid_edge({v, w})); assert(valid_non_edge({v, x})); assert(valid_edge({w, x}));

                                if (callback({y, u, v, w, x})) return true;
                            }
                        }
                    }

                    // Necktie
                    X = neighbors(w) & neighbors(v) & non_neighbors(u);
                    for (Vertex x : Graph::iterate(X)) {
                        Y = neighbors(u) & non_neighbors(v) & non_neighbors(w) & non_neighbors(x);
                        for (Vertex y : Graph::iterate(Y)) {

                            assert(valid_edge({y, u})); assert(valid_non_edge({y, v})); assert(valid_non_edge({y, w})); assert(valid_non_edge({y, x}));
                            assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_non_edge({u, x})); assert(valid_edge({v, w})); assert(valid_edge({v, x})); assert(valid_edge({w, x}));

                            if (callback({y, u, v, w, x})) return true;
                        }
                    }

                    X = neighbors(w) & non_neighbors(v) & non_neighbors(u);
                    for (Vertex x : Graph::iterate(X)) {
                        Y = neighbors(u) & neighbors(v) & non_neighbors(w) & non_neighbors(x);
                        for (Vertex y : Graph::iterate(Y)) {

                            assert(valid_edge({y, u})); assert(valid_edge({y, v})); assert(valid_non_edge({y, w})); assert(valid_non_edge({y, x}));
                            assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_non_edge({u, x})); assert(valid_edge({v, w})); assert(valid_non_edge({v, x})); assert(valid_edge({w, x}));

                            if (callback({y, u, v, w, x})) return true;
                        }
                    }

                    // Bowtie
                    X = neighbors(w) & neighbors(v) & non_neighbors(u);
                    Y = neighbors(u) & neighbors(v) & non_neighbors(w);
                    for (Vertex x : Graph::iterate(X)) {
                        for (Vertex y : Graph::iterate(Y)) {
                            if (!valid_non_edge({x, y})) continue;

                            assert(valid_edge({y, u})); assert(valid_edge({y, v})); assert(valid_non_edge({y, w})); assert(valid_non_edge({y, x}));
                            assert(valid_edge({u, v})); assert(valid_non_edge({u, w})); assert(valid_non_edge({u, x})); assert(valid_edge({v, w})); assert(valid_edge({v, x})); assert(valid_edge({w, x}));

                            if (callback({y, u, v, w, x})) return true;
                        }
                    }
                }
            }
        }
        return false;
    }
}
