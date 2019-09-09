//
// Created by jonas on 02.07.19.
//

#ifndef CONCEPT_NAIVEP3_H
#define CONCEPT_NAIVEP3_H


#include "../interfaces/FinderI.h"


namespace Finder {
    class NaiveP3 : public FinderI {

    public:
        explicit NaiveP3(const Graph &graph_ref) : FinderI(graph_ref) {}

        bool find(SubgraphCallback callback) override;

        bool find(const Graph &forbidden, SubgraphCallback callback) override;

        bool find_near(VertexPair uv, SubgraphCallback callback) override;

        bool find_near(VertexPair uv, const Graph &forbidden, SubgraphCallback callback) override;

        [[nodiscard]] Options::FSG forbidden_subgraphs() const override { return Options::FSG::P3; }

        [[nodiscard]] std::string name() const override { return "NaiveP3"; }

        void to_yaml(YAML::Emitter &out) const override;

    private:
        template <typename H, typename I>
        bool find(const SubgraphCallback &callback, H valid_edge, I valid_non_edge);

        template <typename H, typename I>
        bool find_near(VertexPair uv, const SubgraphCallback &callback, H valid_edge, I valid_non_edge);
    };
}

#endif //CONCEPT_NAIVEP3_H
