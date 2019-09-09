//
// Created by jonas on 02.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINT_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINT_H


#include "detail/EndpointFindImpl.h"


namespace Finder {
    template <size_t length>
    class Endpoint : public FinderI {
        constexpr static bool with_cycles = length > 3;
    public:
        explicit Endpoint(const Graph &graph_ref) : FinderI(graph_ref) {};

        bool find(SubgraphCallback callback) override {
            return detail::EndpointFindImpl<length, with_cycles>::find(graph, callback,
                    neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
        }

        bool find(const Graph& forbidden, SubgraphCallback callback) override {
            return detail::EndpointFindImpl<length, with_cycles>::find(graph, callback,
                    neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
        }

        bool find_near(VertexPair /*uv*/, SubgraphCallback /*callback*/) override {
            assert(false);
            return false;
        }

        bool find_near(VertexPair /*uv*/, const Graph& /*forbidden*/, SubgraphCallback /*callback*/) override  {
            assert(false);
            return false;
        }
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINT_H
