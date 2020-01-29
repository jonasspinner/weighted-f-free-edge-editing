//
// Created by jonas on 25.01.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ILSMWISSOLVERLOWERBOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ILSMWISSOLVERLOWERBOUND_H


#include "../interfaces/LowerBoundI.h"
#include "../../extern/ils_mwis/src/Graph.h"
#include "../Instance.h"

class ILSMWISSolverLowerBound : public LowerBoundI {
private:
    const Graph &m_graph;
    const VertexPairMap <Cost> &m_costs;
    const VertexPairMap<bool> &m_marked;

public:
    ILSMWISSolverLowerBound(const Instance &instance, const VertexPairMap<bool> &marked,
                            std::shared_ptr <FinderI> finder_ref) :
            LowerBoundI(std::move(finder_ref)), m_graph(instance.graph), m_costs(instance.costs), m_marked(marked) {}

    Cost calculate_lower_bound(Cost k) override;

private:
    ils_mwis::Graph build_instance();
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ILSMWISSOLVERLOWERBOUND_H
