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
            bool solved = true;

            finder->find([&](const Subgraph &subgraph) {
                solved = false;

                for (VertexPair uv : subgraph.vertexPairs())
                    if (!m_forbidden[uv])
                        pairs.push_back(uv);

                return true;
            });

            return {pairs, solved};
        }
    };
}


#endif //CONCEPT_FIRSTEDITABLE_H
