//
// Created by jonas on 02.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINT_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINT_H


#include "FinderI.h"


namespace Finder {
    /**
     * A finder class for paths and cycles with a given length.
     * The implementation details can be found in the FindImpl and FindNearImpl classes.
     */
    template <int length, bool with_cycles>
    class Endpoint : public FinderI {
        static_assert(length > 1);

    public:
        bool find(const Graph& graph, SubgraphCallback callback) override;

        bool find(const Graph& graph, const Graph &forbidden, SubgraphCallback callback) override;

        bool find_near(VertexPair uv, const Graph& graph, SubgraphCallback callback) override;

        bool find_near(VertexPair uv, const Graph& graph, const Graph &forbidden, SubgraphCallback callback) override;

        bool find_with_duplicates(const Graph &graph, const SubgraphCallback &callback) override;

        bool find_with_duplicates(const Graph &graph, const Graph &forbidden, const SubgraphCallback &callback) override;

        bool for_all_conversionless_edits(const Subgraph &subgraph, const VertexPairCallback &callback) const override;

        [[nodiscard]] Options::FSG forbidden_subgraphs() const override;

        [[nodiscard]] std::string name() const override;

        void to_yaml(YAML::Emitter &out) const override;

    };

    using EndpointRecC6P6 = Endpoint<6, true>;
    using EndpointRecC5P5 = Endpoint<5, true>;
    using EndpointRecC4P4 = Endpoint<4, true>;
    using EndpointRecP6 = Endpoint<6, false>;
    using EndpointRecP5 = Endpoint<5, false>;
    using EndpointRecP4 = Endpoint<4, false>;
    using EndpointRecP3 = Endpoint<3, false>;
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINT_H
