#ifndef WEIGHTED_F_FREE_EDGE_EDITING_NPS_MWIS_SOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_NPS_MWIS_SOLVER_H


#include <optional>

#include "LowerBoundI.h"
#include "../../extern/nps_mwis/src/Graph.h"
#include "../forbidden_subgraphs/SubgraphC4P4.h"
#include "../editor/EditState.h"


namespace lower_bound {

    template<Options::FSG SetOfForbiddenSubgraphs>
    class NPS_MWIS_Solver : public LowerBoundI {
    private:
        using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
        using Finder = typename Subgraph::Finder;

        const EditState *m_edit_state;

        Finder finder;
    public:
        explicit NPS_MWIS_Solver(const EditState *edit_state) : m_edit_state(edit_state) {
            assert(m_edit_state);
        }

        Cost calculate_lower_bound(Cost k) override;

    private:
        std::optional<nps_mwis::Graph> build_instance();
    };

}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_NPS_MWIS_SOLVER_H
