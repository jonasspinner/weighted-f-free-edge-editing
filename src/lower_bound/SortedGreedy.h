#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SORTEDGREEDY_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SORTEDGREEDY_H


#include "../graph/Graph.h"
#include "../Instance.h"
#include "LowerBoundI.h"
#include "../forbidden_subgraphs/SubgraphC4P4.h"


namespace lower_bound {
    template<Options::FSG SetOfForbiddenSubgraphs>
    class SortedGreedy : public LowerBoundI {
    private:
        using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
        using Finder = typename Subgraph::Finder;

        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;
        VertexPairMap<bool> m_used_in_bound;
        Finder m_finder;
    public:

        SortedGreedy(const Instance &instance, const VertexPairMap<bool> &marked) :
                m_graph(instance.graph), m_costs(instance.costs), m_marked(marked),
                m_used_in_bound(m_graph.size()) {}

        Cost calculate_lower_bound(Cost k) override;

        [[nodiscard]] const auto &used_in_bound() const {
            return m_used_in_bound;
        }

        std::tuple<Cost, VertexPairMap<bool>, std::vector<Subgraph>, std::vector<VertexPair>>
        calculate_lower_bound_and_packing();
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SORTEDGREEDY_H
