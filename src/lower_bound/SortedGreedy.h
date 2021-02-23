#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SORTEDGREEDY_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SORTEDGREEDY_H


#include "LowerBoundI.h"
#include "../forbidden_subgraphs/subgraphs.h"
#include "../editor/EditState.h"


namespace lower_bound {
    template<Options::FSG SetOfForbiddenSubgraphs>
    class SortedGreedy : public LowerBoundI {
    private:
        using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
        using Finder = typename Subgraph::Finder;

        const EditState *m_edit_state;

        VertexPairMap<bool> m_used_in_bound;
        Finder m_finder;
    public:

        explicit SortedGreedy(const EditState *edit_state) :
                m_edit_state(edit_state), m_used_in_bound(m_edit_state->graph().size()) {
            assert(m_edit_state);
        }

        Cost calculate_lower_bound(Cost k) override;

        [[nodiscard]] const auto &used_in_bound() const {
            return m_used_in_bound;
        }

        std::tuple<Cost, VertexPairMap<bool>, std::vector<Subgraph>, std::vector<VertexPair>>
        calculate_lower_bound_and_packing();
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SORTEDGREEDY_H
