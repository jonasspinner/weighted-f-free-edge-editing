#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GREEDY_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GREEDY_H


#include "LowerBoundI.h"
#include "../Instance.h"
#include "../forbidden_subgraphs/SubgraphC4P4.h"


namespace lower_bound {
    template <Options::FSG SetOfForbiddenSubgraphs>
    class Greedy : public LowerBoundI {
    private:
        using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
        using Finder = typename Subgraph::Finder;

        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;
        VertexPairMap<bool> m_used_in_bound;

        Finder finder;
    public:

        Greedy(const Instance &instance, const VertexPairMap<bool> &marked) :
                m_graph(instance.graph), m_costs(instance.costs),
                m_marked(marked), m_used_in_bound(m_graph.size()) {}

        Cost calculate_lower_bound(Cost /*k*/) override;
    };

}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_GREEDY_H
