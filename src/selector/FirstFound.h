#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FIRSTEDITABLE_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FIRSTEDITABLE_H


#include "../graph/VertexPairMap.h"
#include "SelectorI.h"


namespace Selector {

    class FirstFound : public SelectorI {
    private:
        const Graph &m_graph;
        const VertexPairMap<bool> &m_marked;

        std::shared_ptr<FinderI> finder;
    public:
        explicit FirstFound(std::shared_ptr<FinderI> finder_ptr, const Graph &graph, const VertexPairMap<bool> &marked)
                : m_graph(graph), m_marked(marked), finder(std::move(finder_ptr)) {}

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
                    if (!m_marked[uv])
                        problem.pairs.push_back(uv);

                return true;
            });

            return problem;
        }
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FIRSTEDITABLE_H
