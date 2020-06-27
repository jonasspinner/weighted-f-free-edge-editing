#ifndef WEIGHTED_F_FREE_EDGE_EDITING_NPS_MWIS_SOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_NPS_MWIS_SOLVER_H


#include <optional>

#include "LowerBoundI.h"
#include "../../extern/nps_mwis/src/Graph.h"
#include "../Instance.h"
#include "../forbidden_subgraphs/SubgraphC4P4.h"


namespace lower_bound {

    template<Options::FSG SetOfForbiddenSubgraphs>
    class NPS_MWIS_Solver : public LowerBoundI {
    private:
        using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
        using Finder = typename Subgraph::Finder;

        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;

        Finder finder;
    public:
        NPS_MWIS_Solver(const Instance &instance, const VertexPairMap<bool> &marked) :
                m_graph(instance.graph), m_costs(instance.costs), m_marked(marked) {}

        Cost calculate_lower_bound(Cost k) override;

    private:
        std::optional<nps_mwis::Graph> build_instance();
    };

}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_NPS_MWIS_SOLVER_H
