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

    public:
        explicit FirstEditable(std::shared_ptr<FinderI> finder_ptr,
                               const VertexPairMap<bool> &forbidden) : SelectorI(std::move(finder_ptr)),
                                                                       m_forbidden(forbidden) {}

        Problem result(Cost /*k*/) override {
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


#endif //CONCEPT_FIRSTEDITABLE_H
