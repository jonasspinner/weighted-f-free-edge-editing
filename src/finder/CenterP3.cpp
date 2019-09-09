//
// Created by jonas on 23.07.19.
//


#include "CenterP3.h"
#include "../Configuration.h"


namespace Finder {
    /**
     * Calls callback for all P_3's.
     *
     * @param callback
     * @return
     */
    bool CenterP3::find(SubgraphCallback callback) {
        return find(callback, neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
    }

    /**
     * Calls callback for all P_3's. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
     *
     * @param forbidden
     * @param callback
     * @return
     */
    bool CenterP3::find(const Graph &forbidden, SubgraphCallback callback) {
        return find(callback, neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    /**
     * Calls callback for all P_3's having both u and v as vertices.
     *
     * @param uv
     * @param callback
     * @return
     */
    bool CenterP3::find_near(VertexPair uv, SubgraphCallback callback) {
        return find_near(uv, callback, neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
    }

    /**
     * Calls callback for all P_3's having both u and v as vertices. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
     *
     * @param uv
     * @param forbidden
     * @param callback
     * @return
     */
    bool CenterP3::find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) {
        return find_near(uv, callback, neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    void CenterP3::to_yaml(YAML::Emitter &out) const {
        using namespace YAML;
        out << BeginMap;
        out << Key << "name" << Value << "CenterP3";
        out << Key << "forbidden_subgraphs" << Value << Options::FSG::P3;
        out << EndMap;
    }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
    /**
     * Find paths of length 3. Each path is listed exactly once.
     */
    template<typename F, typename G, typename H, typename I>
    bool CenterP3::find(const SubgraphCallback& callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
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
    bool CenterP3::find_near(VertexPair uv, const SubgraphCallback& callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {

        auto [u, v] = uv;

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
}
