//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_FIRSTEDITABLE_H
#define CONCEPT_FIRSTEDITABLE_H


#include "../graph/VertexPairMap.h"
#include "../interfaces/SelectorI.h"


namespace Selector {

    class FirstEditable : public SelectorI {
    private:
        const VertexPairMap<bool> &forbidden;

        class State : public StateI {
            std::unique_ptr<StateI> copy() override {
                return std::make_unique<State>(*this);
            }
        };

    public:
        explicit FirstEditable(const Graph &graph, const VertexPairMap<Cost> &weights, const std::shared_ptr<FinderI> &finder,
                               const VertexPairMap<bool> &forbidden) : SelectorI(finder),
                                                                       forbidden(forbidden) {}

        Problem result(Cost k) override {
            std::vector<VertexPair> pairs;

            bool found = false;
            finder->find([&](const Subgraph &subgraph) {
                subgraph.for_all_unmarked_vertex_pairs(forbidden, [&](VertexPair uv) {
                    pairs.push_back(uv);
                    return false;
                });
                found |= true;
                return !pairs.empty();
            });

            return {pairs, !found};
        }

        std::unique_ptr<StateI> initialize(Cost k) override { return std::make_unique<State>(); }
        void before_mark_and_edit(StateI *state, VertexPair uv) override {}
        void after_mark_and_edit(StateI *state, VertexPair uv) override {}
        void before_mark(StateI *state, VertexPair uv) override {}
        void after_mark(StateI *state, VertexPair uv) override {}
        void before_edit(StateI *state, VertexPair uv) override {}
        void after_edit(StateI *state, VertexPair uv) override {}
        void after_unmark(StateI *state, VertexPair uv) override {}
    };
}


#endif //CONCEPT_FIRSTEDITABLE_H
