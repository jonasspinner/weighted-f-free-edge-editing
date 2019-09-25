//
// Created by jonas on 03.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_CENTERC4P4_H
#define WEIGHTED_F_FREE_EDGE_EDITING_CENTERC4P4_H


#include "../interfaces/FinderI.h"


namespace Finder {
    class CenterC4P4 : public FinderI {

        Graph::AdjRow V;
        Graph::AdjRow A;
        Graph::AdjRow B;

    public:
        explicit CenterC4P4(const Graph &graph_ref) : FinderI(graph_ref), V(graph.size()), A(graph.size()), B(graph.size()) {}

        bool find(SubgraphCallback callback) override;

        bool find(const Graph &forbidden, SubgraphCallback callback) override;

        bool find_near(VertexPair uv, SubgraphCallback callback) override;

        bool find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) override;

        [[nodiscard]] Options::FSG forbidden_subgraphs() const override { return Options::FSG::C4P4; }

        [[nodiscard]] std::string name() const override { return "CenterC4P4"; }

        void to_yaml(YAML::Emitter &out) const override;

    private:
        template<typename F, typename G, typename H, typename I>
        bool find(const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge);

        template<typename F, typename G, typename H, typename I>
        bool find_near(VertexPair uv, const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge,
                       I valid_non_edge);
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_CENTERC4P4_H
