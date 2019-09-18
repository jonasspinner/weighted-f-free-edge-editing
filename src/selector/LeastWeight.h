//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_LEASTWEIGHT_H
#define CONCEPT_LEASTWEIGHT_H

#include "../graph/VertexPairMap.h"

namespace Selector {
    class LeastWeight : public SelectorI {
    private:
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;

    public:
        LeastWeight(const VertexPairMap<Cost> &costs, std::shared_ptr<FinderI> finder_ptr,
                    const VertexPairMap<bool> &marked) : SelectorI(std::move(finder_ptr)), m_costs(costs),
                                                                     m_marked(marked) {}

        Problem result(Cost /*k*/) override {
            Subgraph min_subgraph{};
            Cost min_subgraph_cost = invalid_cost;
            bool solved = true;
            bool unsolveable = false;

            finder->find([&](Subgraph &&subgraph) {
                solved = false;

                Cost subgraph_cost = get_subgraph_cost(subgraph, m_marked, m_costs);

                if (subgraph_cost == invalid_cost) {
                    // if at least one forbidden subgraph has only marked vertex pairs, the problem is not solvable.
                    unsolveable = true;
                } else  if (subgraph_cost < min_subgraph_cost) {
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


#endif //CONCEPT_LEASTWEIGHT_H
