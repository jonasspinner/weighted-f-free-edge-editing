//
// Created by jonas on 29.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GREEDYLOWERBOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GREEDYLOWERBOUND_H


#include "../interfaces/LowerBoundI.h"
#include "../Instance.h"
#include "../graph/Subgraph.h"

namespace LowerBound {
    class GreedyLowerBound : public LowerBoundI {
    private:
        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;
    public:

        GreedyLowerBound(const Instance &instance, const VertexPairMap<bool> &marked,
                         std::shared_ptr<FinderI> finder_ref) : LowerBoundI(std::move(finder_ref)),
                                                                m_graph(instance.graph),
                                                                m_costs(instance.costs), m_marked(marked) {}

        /**
         * Calculates a lower bound on the costs required to solve the current instance.
         *
         * Greedily inserts forbidden subgraphs into the bound. Higher minimum editing costs are preferred. A subgraph
         * is not inserted if it shares an editable vertex pair with a subgraph already in the bound.
         *
         * @param k Not used
         * @return A lower bound on the costs required to solve the current instance.
         */
        Cost result(Cost /*k*/) override {
            /*
            // Find all forbidden subgraphs with editable vertex pairs
            // The cost for a single forbidden subgraph is the minimum edit cost for an editable vertex pair
            std::vector<std::pair<Cost, Subgraph>> subgraphs;
            finder->find([&](Subgraph &&subgraph) {
                Cost min_cost = cost(subgraph, m_marked, m_costs);
                subgraphs.emplace_back(min_cost, std::move(subgraph));
                return false;
            });

            // Sort subgraphs with decreasing costs
            std::sort(subgraphs.begin(), subgraphs.end(),
                      [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

            if (!subgraphs.empty() && subgraphs[0].first == invalid_cost) return 0;

            Cost bound_size = 0;
            VertexPairMap<bool> is_in_bound(m_graph.size(), false);

            // Insert forbidden subgraphs with decreasing minimum edit cost into the bound
            // Only insert a subgraph if it does not share an editable vertex pair with a subgraph already in the bound
            for (const auto&[cost, subgraph] : subgraphs) {
                bool touches_bound = false;
                for (VertexPair uv : subgraph.vertexPairs()) {
                    if (!m_marked[uv] && is_in_bound[uv]) touches_bound = true;
                }

                if (!touches_bound) {
                    bound_size += cost;
                    for (VertexPair uv : subgraph.vertexPairs()) {
                        if (!m_marked[uv])
                        is_in_bound[uv] = true;
                    }
                }
            }
            */

            Cost bound_size = 0;
            VertexPairMap<bool> is_in_bound(m_graph.size(), false);

            finder->find([&](Subgraph &&subgraph) {
                bool touches_bound = false;
                for (VertexPair uv : subgraph.vertexPairs()) {
                    if (!m_marked[uv] && is_in_bound[uv])
                        touches_bound = true;
                }
                if (!touches_bound) {
                    bound_size += get_subgraph_cost(subgraph, m_marked, m_costs);
                    for (VertexPair uv : subgraph.vertexPairs()) {
                        if (!m_marked[uv])
                            is_in_bound[uv] = true;
                    }
                }
                return false;
            });

            return bound_size;
        }
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_GREEDYLOWERBOUND_H
