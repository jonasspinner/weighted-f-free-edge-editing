//
// Created by jonas on 02.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINT_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINT_H


#include "../interfaces/FinderI.h"


namespace Finder {
    template <int length, bool with_cycles>
    class Endpoint : public FinderI {
        static_assert(length > 1);

    public:
        /**
         * A finder class for paths and cycles with a given length.
         * The implementation details can be found in the FindImpl and FindNearImpl classes.
         *
         * @param graph_ref A reference to the graph.
         */
        explicit Endpoint(const Graph &graph_ref) : FinderI(graph_ref) {}

        bool find(SubgraphCallback callback) override;

        bool find(const Graph &forbidden, SubgraphCallback callback) override;

        bool find_near(VertexPair uv, SubgraphCallback callback) override;

        bool find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) override;

        [[nodiscard]] Options::FSG forbidden_subgraphs() const override;

        [[nodiscard]] std::string name() const override;

        void to_yaml(YAML::Emitter &out) const override;

    };

    using EndpointRecC5P5 = Endpoint<5, true>;
    using EndpointRecC4P4 = Endpoint<4, true>;
    using EndpointRecP3 = Endpoint<3, false>;
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ENDPOINT_H
