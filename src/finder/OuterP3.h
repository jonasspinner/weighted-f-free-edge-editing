//
// Created by jonas on 02.11.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_OUTERP3_H
#define WEIGHTED_F_FREE_EDGE_EDITING_OUTERP3_H


#include "../interfaces/FinderI.h"


namespace Finder {
    class OuterP3 : public FinderI {

    public:
        bool find(const Graph& graph, SubgraphCallback callback) override;

        bool find(const Graph& graph, const Graph &forbidden, SubgraphCallback callback) override;

        bool find_near(VertexPair uv, const Graph& graph, SubgraphCallback callback) override;

        bool find_near(VertexPair uv, const Graph& graph, const Graph &forbidden, SubgraphCallback callback) override;

        [[nodiscard]] Options::FSG forbidden_subgraphs() const override { return Options::FSG::P3; }

        [[nodiscard]] std::string name() const override { return "OuterP3"; }

        void to_yaml(YAML::Emitter &out) const override;

    private:
        template <typename F, typename G, typename H, typename I>
        bool find(const Graph& graph, const SubgraphCallback &callback, F neighbor, G non_neighbors, H valid_edge, I valid_non_edge);

        template <typename F, typename G, typename H, typename I>
        bool find_near(VertexPair uv, const Graph& graph, const SubgraphCallback &callback, F neighbor, G non_neighbors, H valid_edge, I valid_non_edge);
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_OUTERP3_H
