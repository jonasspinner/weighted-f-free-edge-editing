//
// Created by jonas on 02.11.19.
//

#include "OuterP3.h"


namespace Finder {
    /**
     * Calls callback for all P_3's.
     *
     * @param graph
     * @param callback
     * @return
     */
    bool OuterP3::find(const Graph& graph, SubgraphCallback callback) {
        return find(graph, callback, neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
    }

    /**
     * Calls callback for all P_3's. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
     *
     * @param graph
     * @param forbidden
     * @param callback
     * @return
     */
    bool OuterP3::find(const Graph& graph, const Graph &forbidden, SubgraphCallback callback) {
        return find(graph, callback, neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    /**
     * Calls callback for all P_3's having both u and v as vertices.
     *
     * @param uv
     * @param graph
     * @param callback
     * @return
     */
    bool OuterP3::find_near(VertexPair uv, const Graph& graph, SubgraphCallback callback) {
        return find_near(uv, graph, callback, neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
    }

    /**
     * Calls callback for all P_3's having both u and v as vertices. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
     *
     * @param uv
     * @param graph
     * @param forbidden
     * @param callback
     * @return
     */
    bool OuterP3::find_near(VertexPair uv, const Graph& graph, const Graph &forbidden, SubgraphCallback callback) {
        return find_near(uv, graph, callback, neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    void OuterP3::to_yaml(YAML::Emitter &out) const {
        using namespace YAML;
        out << BeginMap;
        out << Key << "name" << Value << "OuterP3";
        out << Key << "forbidden_subgraphs" << Value << Options::FSG::P3;
        out << EndMap;
    }

    template<typename F, typename G, typename H, typename I>
    bool OuterP3::find(const Graph& graph, const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
        Graph::AdjRow Y(graph.size()), Z(graph.size());

        /** P_3: <x, y, z> **/
        for (Vertex x : graph.vertices()) {
            // V - N(x) - {x}
            Y = non_neighbors(x);
            for (Vertex y : Graph::iterate(Y)) {
                if (y >= x) continue;

                Z = neighbors(y);
                Z &= neighbors(x);
                for (Vertex z : Graph::iterate(Z)) {
                    assert(y != z); assert(y != x); assert(z != x);

                    assert(valid_edge({y, z})); assert(valid_non_edge({y, x})); assert(valid_edge({z, x}));
                    if (callback(Subgraph{y, z, x})) return true;
                }
            }
        }
        return false;
    }

    template<typename F, typename G, typename H, typename I>
    bool OuterP3::find_near(VertexPair uv, const Graph& graph, const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
        return find(graph, [&](Subgraph &&subgraph){
            const auto &S = subgraph;
            if ((S[0] == uv.u || S[1] == uv.u || S[2] == uv.u) &&
                (S[0] == uv.v || S[1] == uv.v || S[2] == uv.v)) {
                if (callback(std::move(subgraph))) return true;
            }
            return false;
        }, neighbors, non_neighbors, valid_edge, valid_non_edge);
    }

    bool OuterP3::find_with_duplicates(const Graph &graph, const FinderI::SubgraphCallback &callback) {
        return find(graph, callback);
    }

    bool OuterP3::find_with_duplicates(const Graph &graph, const Graph &forbidden,
                                        const FinderI::SubgraphCallback &callback) {
        return find(graph, forbidden, callback);
    }

    bool OuterP3::for_all_conversionless_edits(const Subgraph &subgraph,
                                                const FinderI::VertexPairCallback &callback) const {
        for (auto uv : subgraph.vertexPairs()) {
            if (callback(uv))
                return true;
        }
        return false;
    }
}
