//
// Created by jonas on 29.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_NAIVEC4P4_H
#define WEIGHTED_F_FREE_EDGE_EDITING_NAIVEC4P4_H


#include "FinderI.h"


namespace Finder {
    class NaiveC4P4 : public FinderI {
    public:
        /**
         * Calls callback for all P_4's and C_4's.
         *
         * @param graph
         * @param callback
         * @return
         */
        bool find(const Graph& graph, SubgraphCallback callback) override {
            return find(graph, callback, valid_edge(graph), valid_non_edge(graph));
        }

        /**
         * Calls callback for all P_4's and C_4's. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
         *
         * @param graph
         * @param forbidden
         * @param callback
         * @return
         */
        bool find(const Graph& graph, const Graph &forbidden, SubgraphCallback callback) override {
            return find(graph, callback, valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
        }

        /**
         * Calls callback for all P_4's and C_4's having both u and v as vertices.
         *
         * @param uv
         * @param graph
         * @param callback
         * @return
         */
        bool find_near(VertexPair uv, const Graph& graph, SubgraphCallback callback) override {
            return find_near(uv, graph, callback, valid_edge(graph), valid_non_edge(graph));
        };

        /**
         * Calls callback for all P_4's and C_4's having both u and v as vertices. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
         *
         * @param uv
         * @param graph
         * @param forbidden
         * @param callback
         * @return
         */
        bool find_near(VertexPair uv, const Graph& graph, const Graph &forbidden, SubgraphCallback callback) override {
            return find_near(uv, graph, callback, valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
        }

        [[nodiscard]] Options::FSG forbidden_subgraphs() const override { return Options::FSG::C4P4; }

        [[nodiscard]] std::string name() const override { return "NaiveC4P4"; }

        void to_yaml(YAML::Emitter &out) const override {
            using namespace YAML;
            out << BeginMap;
            out << Key << "name" << Value << "NaiveC4P4";
            out << Key << "forbidden_subgraphs" << Value << Options::FSG::C4P4;
            out << EndMap;
        }

    private:

        template <typename H, typename I>
        bool find(const Graph& graph, const SubgraphCallback& callback, H valid_edge, I valid_non_edge) {

            for (Vertex u : graph.vertices()) {
                for (Vertex v : graph.vertices()) {
                    for (Vertex w : graph.vertices()) {
                        for (Vertex x : graph.vertices()) {
                            if (u != v && u != w && u != x && v != w && v != x && w != x) {
                                /* u-v-w-x */
                                if (valid_edge({u, v}) && valid_non_edge({u, w}) && valid_edge({v, w}) && valid_non_edge({v, x}) && valid_edge({w, x})) {
                                    if (valid_non_edge({u, x})) {
                                        // P_4
                                        if (u < x) // p_1 < p_k
                                            if (callback(Subgraph({u, v, w, x}))) return true;
                                    } else if (valid_edge({u, x})) {
                                        // C_4
                                        if (u < std::min({v, w, x})) // p_1 is smallest
                                            if (v < x) // p_2 < p_k
                                                if (callback(Subgraph({u, v, w, x}))) return true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return false;
        }

        template <typename H, typename I>
        bool find_near(VertexPair uv, const Graph& graph, const SubgraphCallback& callback, H valid_edge, I valid_non_edge) {

            return find(graph, [&](Subgraph &&subgraph) {
                auto vertices = subgraph.vertices();

                bool has_u = std::any_of(vertices.begin(), vertices.end(), [&](Vertex x) { return x == uv.u; });
                bool has_v = std::any_of(vertices.begin(), vertices.end(), [&](Vertex x) { return x == uv.v; });

                if (has_u && has_v) {
                    return callback(std::move(subgraph));
                } else {
                    return false;
                }
            }, valid_edge, valid_non_edge);
        }

    };
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_NAIVEC4P4_H
