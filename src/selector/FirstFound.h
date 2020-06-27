#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FIRSTEDITABLE_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FIRSTEDITABLE_H


#include "../graph/VertexPairMap.h"
#include "SelectorI.h"
#include "../forbidden_subgraphs/SubgraphC4P4.h"


namespace selector {

    template<Options::FSG SetOfForbiddenSubgraphs>
    class FirstFound : public SelectorI {
    private:
        using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
        using Finder = typename Subgraph::Finder;

        const Graph &m_graph;
        const VertexPairMap<bool> &m_marked;

        Finder finder;
    public:
        explicit FirstFound(const Graph &graph, const VertexPairMap<bool> &marked)
                : m_graph(graph), m_marked(marked) {}

        /**
         * Select the first forbidden subgraph found.
         *
         * @return
         */
        Problem select_problem(Cost /*k*/) override {
            Problem problem;
            problem.solved = true;

            finder.find(m_graph, [&](Subgraph subgraph) {
                problem.solved = false;

                for (VertexPair uv : subgraph.vertex_pairs())
                    if (!m_marked[uv])
                        problem.pairs.push_back(uv);

                return true;
            });

            return problem;
        }
    };

    template class FirstFound<Options::FSG::C4P4>;
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FIRSTEDITABLE_H
