//
// Created by jonas on 31.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H


#include "../interfaces/FinderI.h"
#include "center/find.h"
#include "center/find_near.h"


namespace detail {

    template <int length>
    class Center : public FinderI {
        static_assert(length > 1);

    public:
        /**
         * A finder class for paths and cycles with a given length.
         * The implementation details can be found in the FindImpl and FindNearImpl classes.
         *
         * @param graph_ref A reference to the graph.
         */
        explicit Center(const Graph &graph_ref) : FinderI(graph_ref) {}

        /**
         * Calls callback for all paths and cycles of the given length.
         *
         * @param callback
         * @return
         */
        bool find(SubgraphCallback callback) override {
            return detail::FindImpl<length, (length > 3)>::find(graph, callback,
                    neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
        }

        /**
         * Calls callback for all paths and cycles of the given length. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
         *
         * @param forbidden
         * @param callback
         * @return
         */
        bool find(const Graph &forbidden, SubgraphCallback callback) override {
            return detail::FindImpl<length, (length > 3)>::find(graph, callback,
                    neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
        }

        /**
         * Calls callback for all paths and cycles of the given length having both u and v as vertices.
         *
         * @param uv
         * @param callback
         * @return
         */
        bool find_near(VertexPair uv, SubgraphCallback callback) override {
            return detail::FindNearImpl<length, (length > 3)>::find_near(graph, uv, callback,
                    neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
        }

        /**
         * Calls callback for all paths and cycles of the given length having both u and v as vertices. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
         *
         * @param uv
         * @param forbidden
         * @param callback
         * @return
         */
        bool find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) override  {
            return detail::FindNearImpl<length, (length > 3)>::find_near(graph, uv, callback,
                    neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
        }

    };
}

using CenterRecC4P4 = detail::Center<4>;
using CenterRecP3 = detail::Center<3>;

#endif //WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H
