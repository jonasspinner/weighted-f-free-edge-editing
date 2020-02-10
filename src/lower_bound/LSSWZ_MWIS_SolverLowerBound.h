//
// Created by jonas on 06.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LSSWZ_MWIS_SOLVERLOWERBOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LSSWZ_MWIS_SOLVERLOWERBOUND_H

#include "../interfaces/LowerBoundI.h"
#include "../Instance.h"

class LSSWZ_MWIS_SolverLowerBound : public LowerBoundI {
private:
    const Graph &m_graph;
    const VertexPairMap<Cost> &m_costs;
    const VertexPairMap<bool> &m_marked;

public:
    LSSWZ_MWIS_SolverLowerBound(const Instance &instance, const VertexPairMap<bool> &marked,
                                std::shared_ptr<FinderI> finder_ref) :
            LowerBoundI(std::move(finder_ref)), m_graph(instance.graph), m_costs(instance.costs), m_marked(marked) {}

    Cost calculate_lower_bound(Cost k) override;

private:
    static std::pair<Graph, std::vector<Cost>>
    build_instance(FinderI &finder, const Graph &graph, const VertexPairMap<bool> &marked,
                   const VertexPairMap<Cost> &costs);
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_LSSWZ_MWIS_SOLVERLOWERBOUND_H
