//
// Created by jonas on 29.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GREEDYLOWERBOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GREEDYLOWERBOUND_H


#include "LowerBoundI.h"
#include "../Instance.h"


namespace lower_bound {
    class GreedyLowerBound : public LowerBoundI {
    private:
        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;
        VertexPairMap<bool> m_used_in_bound;
    public:

        GreedyLowerBound(const Instance &instance, const VertexPairMap<bool> &marked,
                         std::shared_ptr<FinderI> finder_ref) :
                LowerBoundI(std::move(finder_ref)), m_graph(instance.graph), m_costs(instance.costs),
                m_marked(marked), m_used_in_bound(m_graph.size()) {}

        Cost calculate_lower_bound(Cost /*k*/) override;
    };
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_GREEDYLOWERBOUND_H
