//
// Created by jonas on 25.01.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_NPS_MWIS_SOLVERLOWERBOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_NPS_MWIS_SOLVERLOWERBOUND_H


#include <optional>

#include "../interfaces/LowerBoundI.h"
#include "../../extern/nps_mwis/src/Graph.h"
#include "../Instance.h"

class NPS_MWIS_SolverLowerBound : public LowerBoundI {
private:
    const Graph &m_graph;
    const VertexPairMap <Cost> &m_costs;
    const VertexPairMap<bool> &m_marked;

public:
    NPS_MWIS_SolverLowerBound(const Instance &instance, const VertexPairMap<bool> &marked,
                              std::shared_ptr <FinderI> finder_ref) :
            LowerBoundI(std::move(finder_ref)), m_graph(instance.graph), m_costs(instance.costs), m_marked(marked) {}

    Cost calculate_lower_bound(Cost k) override;

private:
    std::optional<nps_mwis::Graph> build_instance();
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_NPS_MWIS_SOLVERLOWERBOUND_H
