//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_LEASTWEIGHT_H
#define CONCEPT_LEASTWEIGHT_H

#include "../graph/VertexPairMap.h"

namespace Selector {
    class LeastWeight : public SelectorI {
    private:
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_forbidden;

    public:
        explicit LeastWeight(const VertexPairMap<Cost> &costs, std::shared_ptr<FinderI> finder_ptr,
                             const VertexPairMap<bool> &forbidden) : SelectorI(std::move(finder_ptr)), m_costs(costs),
                                                                     m_forbidden(forbidden) {}

        Problem result(Cost /*k*/) override {
            Subgraph min_subgraph{};
            Cost min_subgraph_cost = invalid_cost;


            bool solved = true;

            finder->find([&](Subgraph &&subgraph) {
                solved = false;

                Cost subgraph_cost = cost(subgraph, m_forbidden, m_costs);
                if (subgraph_cost < min_subgraph_cost) {
                    min_subgraph_cost = subgraph_cost;
                    min_subgraph = std::move(subgraph);
                }
                return false;
            });


            std::vector<VertexPair> pairs;

            for (VertexPair uv : min_subgraph.vertexPairs()) {
                if (!m_forbidden[uv])
                    pairs.push_back(uv);
            }

            std::sort(pairs.begin(), pairs.end(),
                      [&](VertexPair uv, VertexPair xy) { return m_costs[uv] < m_costs[xy]; });

            return {pairs, solved};
        }

        void push(Cost /*k*/) override {}

        void pop() override {}

        void before_mark_and_edit(VertexPair) override {}

        void after_mark_and_edit(VertexPair) override {}

        void before_mark(VertexPair) override {}

        void after_mark(VertexPair) override {}

        void before_edit(VertexPair) override {}

        void after_edit(VertexPair) override {}

        void after_unmark(VertexPair) override {}
    };
}


#endif //CONCEPT_LEASTWEIGHT_H
