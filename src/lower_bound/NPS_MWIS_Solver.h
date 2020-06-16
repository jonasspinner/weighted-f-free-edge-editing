#ifndef WEIGHTED_F_FREE_EDGE_EDITING_NPS_MWIS_SOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_NPS_MWIS_SOLVER_H


#include <optional>

#include "LowerBoundI.h"
#include "../../extern/nps_mwis/src/Graph.h"
#include "../Instance.h"


namespace lower_bound {

    class NPS_MWIS_Solver : public LowerBoundI {
    private:
        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;

        std::shared_ptr<FinderI> finder;
    public:
        NPS_MWIS_Solver(const Instance &instance, const VertexPairMap<bool> &marked,
                        std::shared_ptr<FinderI> finder_ref) :
                m_graph(instance.graph), m_costs(instance.costs),
                m_marked(marked), finder(std::move(finder_ref)) {}

        Cost calculate_lower_bound(Cost k) override;

    private:
        std::optional<nps_mwis::Graph> build_instance();
    };

}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_NPS_MWIS_SOLVER_H
