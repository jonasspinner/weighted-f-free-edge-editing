#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GREEDY_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GREEDY_H


#include "LowerBoundI.h"
#include "../forbidden_subgraphs/SubgraphC4P4.h"
#include "../editor/EditState.h"


namespace lower_bound {
    template<Options::FSG SetOfForbiddenSubgraphs>
    class Greedy : public LowerBoundI {
    private:
        using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
        using Finder = typename Subgraph::Finder;

        const EditState *m_edit_state;

        VertexPairMap<bool> m_used_in_bound;

        Finder finder;
    public:

        explicit Greedy(const EditState *edit_state) :
                m_edit_state(edit_state), m_used_in_bound(m_edit_state->graph().size()) {}

        Cost calculate_lower_bound(Cost /*k*/) override;
    };

}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_GREEDY_H
