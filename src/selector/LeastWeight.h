#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LEASTWEIGHT_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LEASTWEIGHT_H

#include "../graph/VertexPairMap.h"
#include "../forbidden_subgraphs/SubgraphC4P4.h"


namespace selector {

    template<Options::FSG SetOfForbiddenSubgraphs>
    class LeastWeight : public SelectorI {
    private:
        using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
        using Finder = typename Subgraph::Finder;

        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;

        Finder finder;
    public:
        LeastWeight(const Graph &graph, const VertexPairMap<Cost> &costs,
                    const VertexPairMap<bool> &marked) :
                m_graph(graph), m_costs(costs), m_marked(marked) {}

        Problem select_problem(Cost /*k*/) override {
            std::optional<Subgraph> min_subgraph;
            Cost min_subgraph_cost = invalid_cost;
            bool unsolveable = false;

            finder.find(m_graph, [&](Subgraph subgraph) {
                Cost subgraph_cost = subgraph.calculate_min_cost(m_costs, m_marked);

                if (subgraph_cost == invalid_cost) {
                    // if at least one forbidden subgraph has only marked vertex pairs, the problem is not solvable.
                    unsolveable = true;
                } else if (subgraph_cost < min_subgraph_cost) {
                    // update subgraph with minimum edit cost
                    min_subgraph_cost = subgraph_cost;
                    min_subgraph = std::move(subgraph);
                }
                return subgraph_iterators::break_if(unsolveable);
            });

            if (unsolveable)
                return {{}, false};
            if (!min_subgraph)
                return {{}, true};

            std::vector<VertexPair> pairs;
            for (VertexPair uv : min_subgraph->vertex_pairs())
                if (!m_marked[uv])
                    pairs.push_back(uv);

            std::sort(pairs.begin(), pairs.end(),
                      [&](VertexPair uv, VertexPair xy) { return m_costs[uv] < m_costs[xy]; });

            return {pairs, false};
        }
    };

    template
    class LeastWeight<Options::FSG::C4P4>;
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_LEASTWEIGHT_H
