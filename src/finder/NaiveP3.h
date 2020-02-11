//
// Created by jonas on 02.07.19.
//

#ifndef CONCEPT_NAIVEP3_H
#define CONCEPT_NAIVEP3_H


#include "../interfaces/FinderI.h"


namespace Finder {
    class NaiveP3 : public FinderI {

    public:
        bool find(const Graph& graph, SubgraphCallback callback) override;

        bool find(const Graph& graph, const Graph &forbidden, SubgraphCallback callback) override;

        bool find_near(VertexPair uv, const Graph& graph, SubgraphCallback callback) override;

        bool find_near(VertexPair uv, const Graph& graph, const Graph &forbidden, SubgraphCallback callback) override;

        [[nodiscard]] Options::FSG forbidden_subgraphs() const override { return Options::FSG::P3; }

        [[nodiscard]] std::string name() const override { return "NaiveP3"; }

        void to_yaml(YAML::Emitter &out) const override;

    private:
        template <typename H, typename I>
        bool find(const Graph& graph, const SubgraphCallback &callback, H valid_edge, I valid_non_edge);

        template <typename H, typename I>
        bool find_near(VertexPair uv, const Graph& graph, const SubgraphCallback &callback, H valid_edge, I valid_non_edge);
    };
}

#endif //CONCEPT_NAIVEP3_H
