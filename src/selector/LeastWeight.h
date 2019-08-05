//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_LEASTWEIGHT_H
#define CONCEPT_LEASTWEIGHT_H

#include "../graph/VertexPairMap.h"

namespace Selector {
    class LeastWeight : public SelectorI {
    private:
        const VertexPairMap<Cost> &weights;
        const VertexPairMap<bool> &forbidden;

        class State : public StateI {
            std::unique_ptr<StateI> copy() override {
                return std::make_unique<State>(*this);
            }
        };

    public:
        explicit LeastWeight(const VertexPairMap<Cost> &weights, const std::shared_ptr<FinderI> &finder,
                             const VertexPairMap<bool> &forbidden) : SelectorI(finder), weights(weights),
                                                                     forbidden(forbidden) {}

        Problem result(StateI &, Cost k) override {
            Subgraph min_subgraph{};
            Cost min_subgraph_cost = std::numeric_limits<Cost>::max();

            bool found = false;
            this->finder->find([&](auto subgraph) {
                found |= true;
                Cost subgraph_cost = std::numeric_limits<Cost>::max();
                subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
                    subgraph_cost = std::min(subgraph_cost, weights[uv]);
                    return false;
                });
                if (subgraph_cost < min_subgraph_cost) {
                    min_subgraph_cost = subgraph_cost;
                    min_subgraph = subgraph;
                }
                return false;
            });

            std::vector<VertexPair> pairs;
            min_subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
                pairs.push_back(uv);
                return false;
            });

            std::sort(pairs.begin(), pairs.end(),
                      [&](const VertexPair &uv, const VertexPair &xy) { return weights[uv] > weights[xy]; });

            return {pairs, !found};
        }

        std::unique_ptr<StateI> initialize(Cost k) override { return std::make_unique<State>(); }

        void before_mark_and_edit(StateI &state, VertexPair uv) override {}

        void after_mark_and_edit(StateI &state, VertexPair uv) override {}

        void before_mark(StateI &state, VertexPair uv) override {}

        void after_mark(StateI &state, VertexPair uv) override {}

        void before_edit(StateI &state, VertexPair uv) override {}

        void after_edit(StateI &state, VertexPair uv) override {}

        void after_unmark(StateI &state, VertexPair uv) override {}
    };
}


#endif //CONCEPT_LEASTWEIGHT_H
