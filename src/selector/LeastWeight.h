#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LEASTWEIGHT_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LEASTWEIGHT_H

#include "../graph/VertexPairMap.h"

namespace Selector {
    class LeastWeight : public SelectorI {
    private:
        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;

        std::shared_ptr<FinderI> finder;
    public:
        LeastWeight(const Graph &graph, const VertexPairMap<Cost> &costs, std::shared_ptr<FinderI> finder_ptr,
                    const VertexPairMap<bool> &marked) :
                    m_graph(graph), m_costs(costs), m_marked(marked), finder(std::move(finder_ptr)) {}

        Problem select_problem(Cost /*k*/) override {
            Subgraph min_subgraph{};
            Cost min_subgraph_cost = invalid_cost;
            bool solved = true;
            bool unsolveable = false;

            finder->find(m_graph, [&](Subgraph &&subgraph) {
                solved = false;

                Cost subgraph_cost = finder->calculate_min_cost(subgraph, m_marked, m_costs);

                if (subgraph_cost == invalid_cost) {
                    // if at least one forbidden subgraph has only marked vertex pairs, the problem is not solvable.
                    unsolveable = true;
                } else if (subgraph_cost < min_subgraph_cost) {
                    // update subgraph with minimum edit cost
                    min_subgraph_cost = subgraph_cost;
                    min_subgraph = std::move(subgraph);
                }
                return unsolveable;
            });

            if (unsolveable)
                return {{}, false};

            std::vector<VertexPair> pairs;
            for (VertexPair uv : min_subgraph.vertexPairs())
                if (!m_marked[uv])
                    pairs.push_back(uv);

            std::sort(pairs.begin(), pairs.end(),
                      [&](VertexPair uv, VertexPair xy) { return m_costs[uv] < m_costs[xy]; });

            return {pairs, solved};
        }
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_LEASTWEIGHT_H
