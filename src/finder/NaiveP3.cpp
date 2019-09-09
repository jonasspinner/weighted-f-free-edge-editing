//
// Created by jonas on 02.07.19.
//


#include "NaiveP3.h"


namespace Finder {
    /**
     * Calls callback for all P_3's.
     *
     * @param callback
     * @return
     */
    bool NaiveP3::find(SubgraphCallback callback) {
        return find(callback, valid_edge(graph), valid_non_edge(graph));
    }

    /**
     * Calls callback for all P_3's. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
     *
     * @param forbidden
     * @param callback
     * @return
     */
    bool NaiveP3::find(const Graph &forbidden, SubgraphCallback callback) {
        return find(callback, valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    /**
     * Calls callback for all P_3's having both u and v as vertices.
     *
     * @param uv
     * @param callback
     * @return
     */
    bool NaiveP3::find_near(VertexPair uv, SubgraphCallback callback) {
        return find_near(uv, callback, valid_edge(graph), valid_non_edge(graph));
    }

    /**
     * Calls callback for all P_3's having both u and v as vertices. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
     *
     * @param uv
     * @param forbidden
     * @param callback
     * @return
     */
    bool NaiveP3::find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) {
        return find_near(uv, callback, valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    void NaiveP3::to_yaml(YAML::Emitter &out) const {
        using namespace YAML;
        out << BeginMap;
        out << Key << "name" << Value << "NaiveP3";
        out << Key << "forbidden_subgraphs" << Value << Options::FSG::P3;
        out << EndMap;
    }

    template<typename H, typename I>
    bool NaiveP3::find(const SubgraphCallback &callback, H valid_edge, I valid_non_edge) {
        for (Vertex u : graph.vertices()) {
            for (Vertex v : graph.vertices()) {
                for (Vertex w : graph.vertices()) {
                    if (u != v && u != w && v != w) { // u, v, w are three distinct vertices
                        if (valid_edge({u, v}) && valid_non_edge({u, w}) && valid_edge({u, w})) { // uvw is a P3
                            if (u < w) // break symmetry and list {u, v, w} only once
                                if (callback(Subgraph{u, v, w})) return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    template<typename H, typename I>
    bool NaiveP3::find_near(VertexPair uv, const SubgraphCallback &callback, H valid_edge, I valid_non_edge) {
        return find([&](Subgraph &&subgraph){
            const auto &S = subgraph;
            if ((S[0] == uv.u || S[1] == uv.u || S[2] == uv.u) &&
                (S[0] == uv.v || S[1] == uv.v || S[2] == uv.v)) {
                if (callback(std::move(subgraph))) return true;
            }
            return false;
        }, valid_edge, valid_non_edge);
    }

}
