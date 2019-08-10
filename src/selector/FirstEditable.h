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
        const VertexPairMap<bool> &m_forbidden;

        class State : public StateI {
            std::unique_ptr<StateI> copy() override {
                return std::make_unique<State>(*this);
            }
        };

    public:
        explicit FirstEditable(std::shared_ptr<FinderI> finder_ptr,
                               const VertexPairMap<bool> &forbidden) : SelectorI(finder_ptr),
                                                                       m_forbidden(forbidden) {}

        Problem result(StateI &, Cost /*k*/) override {
            std::vector<VertexPair> pairs;

            bool found = false;
            finder->find([&](const Subgraph &subgraph) {
                for (VertexPair uv : subgraph.vertexPairs()) {
                    if (!m_forbidden[uv]) {
                        pairs.push_back(uv);
                    }
                }
                found |= true;
                return !pairs.empty();
            });

            return {pairs, !found};
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


#endif //CONCEPT_FIRSTEDITABLE_H
