//
// Created by jonas on 02.11.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_NAIVE_H
#define WEIGHTED_F_FREE_EDGE_EDITING_NAIVE_H


#include "FinderI.h"


namespace Finder {
    /**
     * A finder class for paths and cycles with a given length.
     */
    template <int length, bool with_cycles>
    class Naive : public FinderI {
        static_assert(length > 1);

    public:
        bool find(const Graph& graph, SubgraphCallback callback) override;

        bool find(const Graph& graph, const Graph &forbidden, SubgraphCallback callback) override;

        bool find_near(VertexPair uv, const Graph& graph, SubgraphCallback callback) override;

        bool find_near(VertexPair uv, const Graph& graph, const Graph &forbidden, SubgraphCallback callback) override;

        [[nodiscard]] Options::FSG forbidden_subgraphs() const override;

        [[nodiscard]] std::string name() const override;

        void to_yaml(YAML::Emitter &out) const override;

    };

    using NaiveRecC6P6 = Naive<6, true>;
    using NaiveRecC5P5 = Naive<5, true>;
    using NaiveRecC4P4 = Naive<4, true>;
    using NaiveRecP6 = Naive<6, false>;
    using NaiveRecP5 = Naive<5, false>;
    using NaiveRecP4 = Naive<4, false>;
    using NaiveRecP3 = Naive<3, false>;
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_NAIVE_H
