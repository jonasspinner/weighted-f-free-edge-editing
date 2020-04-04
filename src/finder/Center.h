//
// Created by jonas on 31.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H


#include "FinderI.h"


namespace Finder {

    /**
     * A finder class for paths and cycles with a given length.
     * The implementation details can be found in the FindImpl and FindNearImpl classes.
     */
    template <int length, bool with_cycles>
    class Center : public FinderI {
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

    using CenterRecC6P6 = Center<6, true>;
    using CenterRecC5P5 = Center<5, true>;
    using CenterRecC4P4 = Center<4, true>;
    using CenterRecP6 = Center<6, false>;
    using CenterRecP5 = Center<5, false>;
    using CenterRecP4 = Center<4, false>;
    using CenterRecP3 = Center<3, false>;
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_CENTER_H
