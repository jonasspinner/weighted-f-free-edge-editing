#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GREEDY_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GREEDY_H


#include "LowerBoundI.h"
#include "../Instance.h"


namespace lower_bound {

    class Greedy : public LowerBoundI {
    private:
        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;
        VertexPairMap<bool> m_used_in_bound;

        std::shared_ptr<FinderI> finder;
    public:

        Greedy(const Instance &instance, const VertexPairMap<bool> &marked,
               std::shared_ptr<FinderI> finder_ref) :
                m_graph(instance.graph), m_costs(instance.costs),
                m_marked(marked), m_used_in_bound(m_graph.size()), finder(std::move(finder_ref)) {}

        Cost calculate_lower_bound(Cost /*k*/) override;
    };

}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_GREEDY_H
