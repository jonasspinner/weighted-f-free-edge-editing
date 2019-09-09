//
// Created by jonas on 29.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LOCALGREEDYLOWERBOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LOCALGREEDYLOWERBOUND_H


#include "../interfaces/LowerBoundI.h"
#include "../Instance.h"


namespace LowerBound {
    class LocalGreedyLowerBound : public LowerBoundI {
    private:
        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;
    public:

        LocalGreedyLowerBound(const Instance &instance, const VertexPairMap<bool> &marked,
                              std::shared_ptr<FinderI> finder_ref) : LowerBoundI(std::move(finder_ref)),
                                                                     m_graph(instance.graph),
                                                                     m_costs(instance.costs), m_marked(marked) {}

        Cost result(Cost /*k*/) override;
    };
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_LOCALGREEDYLOWERBOUND_H
