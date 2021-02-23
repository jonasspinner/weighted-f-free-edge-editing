#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FIRSTEDITABLE_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FIRSTEDITABLE_H


#include "../graph/VertexPairMap.h"
#include "SelectorI.h"
#include "../forbidden_subgraphs/subgraphs.h"


namespace selector {

    template<Options::FSG SetOfForbiddenSubgraphs>
    class FirstFound : public SelectorI {
    private:
        using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
        using Finder = typename Subgraph::Finder;

        const EditState *m_edit_state;

        Finder finder;
    public:
        explicit FirstFound(const EditState *edit_state) : m_edit_state(edit_state) {
            assert(edit_state);
        }

        /**
         * Select the first forbidden subgraph found.
         *
         * @return
         */
        Problem select_problem(Cost /*k*/) override {
            Problem problem;
            problem.solved = true;

            finder.find(m_edit_state->graph(), [&](Subgraph subgraph) {
                problem.solved = false;

                for (VertexPair uv : subgraph.vertex_pairs())
                    if (!m_edit_state->is_marked(uv))
                        problem.pairs.push_back(uv);

                return subgraph_iterators::IterationControl::Break;
            });

            return problem;
        }
    };

    template class FirstFound<Options::FSG::C4P4>;
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FIRSTEDITABLE_H
