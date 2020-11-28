#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LSSWZ_MWIS_SOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LSSWZ_MWIS_SOLVER_H


#include "LowerBoundI.h"
#include "../options.h"
#include "../Configuration.h"
#include "../forbidden_subgraphs/SubgraphC4P4.h"
#include "../editor/EditState.h"


namespace lower_bound {

    template<Options::FSG SetOfForbiddenSubgraphs>
    class LSSWZ_MWIS_Solver : public LowerBoundI {
    private:
        using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
        using Finder = typename Subgraph::Finder;

        const EditState *m_edit_state;

        int m_verbosity;

        double time_limit = 10;
        bool disable_reduction = false;

        Finder finder;
    public:
        LSSWZ_MWIS_Solver(const EditState *edit_state, int verbosity) :
                m_edit_state(edit_state), m_verbosity(verbosity) {}

        Cost calculate_lower_bound(Cost k) override;

    private:
        static std::optional<std::pair<Graph, std::vector<Cost>>>
        build_instance(Finder &finder, const Graph &graph, const VertexPairMap<bool> &marked,
                       const VertexPairMap<Cost> &costs);
    };

}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_LSSWZ_MWIS_SOLVER_H
