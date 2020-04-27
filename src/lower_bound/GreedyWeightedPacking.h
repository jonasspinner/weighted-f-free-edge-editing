#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GREEDYWEIGHTEDPACKING_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GREEDYWEIGHTEDPACKING_H


#include <queue>

#include "../definitions.h"
#include "LowerBoundI.h"
#include "../Instance.h"


namespace lower_bound {
    class GreedyWeightedPacking : public LowerBoundI {
    private:
        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;
        VertexPairMap<Cost> m_costs_remaining;
    public:
        GreedyWeightedPacking(const Instance &instance, const VertexPairMap<bool> &marked,
                              std::shared_ptr<FinderI> finder_ref) :
                LowerBoundI(std::move(finder_ref)), m_graph(instance.graph), m_costs(instance.costs), m_marked(marked),
                m_costs_remaining(0) {}

        Cost calculate_lower_bound(Cost k) override {
            // The subgraphs are stored in a vector. The priority queue stores the subgraph costs and an index into the
            // vector.
            std::vector<Subgraph> subgraphs;
            std::priority_queue<std::pair<Cost, size_t>> Q;
            Cost lower_bound = 0;
            Cost max_min_cost = std::numeric_limits<Cost>::min();

            // Initialize remaining costs with the editing costs.
            m_costs_remaining = m_costs;


            bool early_exit = finder->find_with_duplicates(m_graph, [&](Subgraph &&subgraph) {
                Cost initial_min_cost = finder->calculate_min_cost(subgraph, m_marked, m_costs_remaining);

                subgraphs.emplace_back(std::move(subgraph));
                Q.emplace(initial_min_cost, subgraphs.size() - 1);

                max_min_cost = std::max(max_min_cost, initial_min_cost);
                return max_min_cost > k;
            });

            // If a single subgraph has an editing cost larger than k, or has an invalid edititing cost (i.e. only has
            // edits that are either marked, or lead to a conversion to another forbidden subgraph), the current
            // instance is no longer solvable.
            if (early_exit)
                return max_min_cost;


            while (!Q.empty() && lower_bound < k) {
                auto[cost, index] = Q.top(); Q.pop();

                // m_cost_remaining may be updated after cost has been calculated.
                Cost current_min_cost = finder->calculate_min_cost(subgraphs[index], m_marked, m_costs_remaining);

                if (current_min_cost == 0) {
                    // The subgraph cannot be inserted.
                    continue;
                } else if (current_min_cost == cost) {
                    // The editing cost is maximal from all remaining subgraphs in the queue. Increase the lower bound
                    // and update the remaining cost matrix.
                    lower_bound += current_min_cost;

                    finder->for_all_conversionless_edits(subgraphs[index], [&](VertexPair uv) {
                        m_costs_remaining[uv] -= current_min_cost;
                        return false;
                    });
                } else {
                    // The editing cost is no longer up to date. The subgraph may be inserted in the future.
                    Q.emplace(current_min_cost, index);
                }
            }

            return lower_bound;
        }
    };
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_GREEDYWEIGHTEDPACKING_H
