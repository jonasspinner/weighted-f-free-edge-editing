//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_FIRSTEDITABLE_H
#define CONCEPT_FIRSTEDITABLE_H


#include "../graph/VertexPairMap.h"
#include "SelectorI.h"


namespace Selector {

    class FirstFound : public SelectorI {
    private:
        const Graph &m_graph;
        const VertexPairMap<bool> &m_forbidden;

    public:
        explicit FirstFound(std::shared_ptr<FinderI> finder_ptr, const Graph &graph,
                            const VertexPairMap<bool> &forbidden) : SelectorI(std::move(finder_ptr)), m_graph(graph),
                                                                    m_forbidden(forbidden) {}

        /**
         * Select the first forbidden subgraph found.
         *
         * @return
         */
        Problem select_problem(Cost /*k*/) override {
            Problem problem;
            problem.solved = true;

            finder->find(m_graph, [&](const Subgraph &subgraph) {
                problem.solved = false;

                for (VertexPair uv : subgraph.vertexPairs())
                    if (!m_forbidden[uv])
                        problem.pairs.push_back(uv);

                return true;
            });

            return problem;
        }
    };
}


#endif //CONCEPT_FIRSTEDITABLE_H
