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

        class State : public StateI {
            std::unique_ptr<StateI> copy() override {
                return std::make_unique<State>(*this);
            }
        };

    public:
        explicit LeastWeight(const VertexPairMap<Cost> &costs, std::shared_ptr<FinderI> finder_ptr,
                             const VertexPairMap<bool> &forbidden) : SelectorI(std::move(finder_ptr)), m_costs(costs),
                                                                     m_forbidden(forbidden) {}

        Problem result(StateI &, Cost /*k*/) override {
            Subgraph min_subgraph{};
            Cost min_subgraph_cost = invalid_cost;

            bool no_subgraphs = true;
            this->finder->find([&](Subgraph &&subgraph) {
                no_subgraphs = false;
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

            return {pairs, no_subgraphs};
        }

        std::unique_ptr<StateI> initialize(Cost /*k*/) override { return std::make_unique<State>(); }

        void before_mark_and_edit(StateI &/*state*/, VertexPair /*uv*/) override {}

        void after_mark_and_edit(StateI &/*state*/, VertexPair /*uv*/) override {}

        void before_mark(StateI &/*state*/, VertexPair /*uv*/) override {}

        void after_mark(StateI &/*state*/, VertexPair /*uv*/) override {}

        void before_edit(StateI &/*state*/, VertexPair /*uv*/) override {}

        void after_edit(StateI &/*state*/, VertexPair /*uv*/) override {}

        void after_unmark(StateI &/*state*/, VertexPair /*uv*/) override {}
    };
}


#endif //CONCEPT_LEASTWEIGHT_H
