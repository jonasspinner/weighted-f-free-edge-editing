//
// Created by jonas on 02.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SPLITGRAPH_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SPLITGRAPH_H


#include "../interfaces/FinderI.h"


namespace Finder {
    class SplitGraph : public FinderI {
    private:
        Graph::AdjRow V;
        Graph::AdjRow W;
        Graph::AdjRow X;
        Graph::AdjRow Y;
    public:
        explicit SplitGraph(const Graph &graph_ref) : FinderI(graph_ref), V(graph.size()), W(graph.size()), X(graph.size()), Y(graph.size()) {}

        bool find(SubgraphCallback callback) override;

        bool find(const Graph& forbidden, SubgraphCallback callback) override;

        bool find_near(VertexPair /*uv*/, SubgraphCallback /*callback*/) override;

        bool find_near(VertexPair /*uv*/, const Graph& /*forbidden*/, SubgraphCallback /*callback*/) override;

        [[nodiscard]] Options::FSG forbidden_subgraphs() const override { return Options::FSG::C4_C5_2K2; }

        [[nodiscard]] std::string name() const override { return "SplitGraph"; }

        void to_yaml(YAML::Emitter &out) const override;

    private:
        template <typename F, typename G, typename H, typename I>
        bool find(const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge);

        template <typename F, typename G, typename H, typename I>
        bool find_near(VertexPair uv, const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge);
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SPLITGRAPH_H
