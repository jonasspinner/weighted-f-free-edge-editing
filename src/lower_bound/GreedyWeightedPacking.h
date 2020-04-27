//
// Created by jonas on 08.04.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GREEDYWEIGHTEDPACKING_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GREEDYWEIGHTEDPACKING_H


#include <queue>

class WeightedPacking {
    VertexPairMap<Cost> m_potential;
    Graph m_saturated_vertex_pairs;
    Cost m_total_cost = 0;
    std::unique_ptr<FinderI> m_finder;
    Graph m_graph;
    std::vector<std::pair<Subgraph, Cost>> m_subgraphs_in_packing;

    void insert_subgraph(const Subgraph &subgraph) {
        Cost cost = calculate_min_cost(subgraph);

        m_finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
            assert(!m_saturated_vertex_pairs.hasEdge(uv));
            m_potential[uv] -= cost;
            if (m_potential[uv] == 0) {
                m_saturated_vertex_pairs.setEdge(uv);
            }
            return false;
        });

        m_total_cost += cost;

        m_subgraphs_in_packing.emplace_back(subgraph, cost);
    }

    void remove_subgraph(const std::pair<Subgraph, Cost> &element) {
        const Subgraph &subgraph = element.first;
        Cost cost = element.second;

        m_finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
            m_potential[uv] += cost;
            if (m_potential[uv] > 0) {
                m_saturated_vertex_pairs.clearEdge(uv);
            }
            return false;
        });
    }

    std::vector<Subgraph> neighbors(const Subgraph &subgraph) {
        std::vector<Subgraph> result;

        m_finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
            assert(!m_saturated_vertex_pairs.hasEdge(uv));
            m_finder->find_near_with_duplicates(uv, m_graph, m_saturated_vertex_pairs, [&](Subgraph &&neighbor) {
                result.push_back(std::move(neighbor));
                return false;
            });
            m_saturated_vertex_pairs.setEdge(uv);
            return false;
        });

        m_finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
            assert(m_saturated_vertex_pairs.hasEdge(uv));
            m_saturated_vertex_pairs.clearEdge(uv);
            return false;
        });
        return result;
    }

    void find_one_two_improvement(size_t i) {
        auto old_element = m_subgraphs_in_packing[i];

        remove_subgraph(old_element);

        auto n = neighbors(old_element.first);

        for (size_t a_index = 0; a_index < n.size(); ++a_index) {
            for (size_t b_index = a_index + 1; b_index < n.size(); ++b_index) {

            }
        }
    }

    Cost calculate_min_cost(const Subgraph &subgraph) {
        Cost cost = std::numeric_limits<Cost>::max();
        m_finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
            assert(!m_saturated_vertex_pairs.hasEdge(uv));
            cost = std::min(cost, m_potential[uv]);
            return false;
        });
        return cost;
    }

    bool is_maximal() {
        return m_finder->find(m_graph, m_saturated_vertex_pairs, [](auto) {
            return true;
        });
    }
};

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
                m_costs_remaining(m_costs.size()) {}

        Cost calculate_lower_bound(Cost k) override {
            std::priority_queue<std::pair<Cost, Subgraph>> Q;
            Cost lower_bound = 0;
            Cost max_min_cost = std::numeric_limits<Cost>::min();

            m_costs_remaining = m_costs;

            bool early_exit = finder->find_with_duplicates(m_graph, [&](Subgraph &&subgraph) {
                Cost initial_min_cost = finder->calculate_min_cost(subgraph, m_marked, m_costs_remaining);
                Q.emplace(initial_min_cost, std::move(subgraph));
                max_min_cost = std::max(max_min_cost, initial_min_cost);
                return max_min_cost > k;
            });

            if (early_exit)
                return max_min_cost;


            while (!Q.empty() && lower_bound < k) {
                auto x = Q.top();
                Q.pop();

                Cost current_min_cost = finder->calculate_min_cost(x.second, m_marked, m_costs_remaining);

                if (current_min_cost == 0) {
                    continue;
                } else if (current_min_cost == x.first) {
                    lower_bound += current_min_cost;

                    finder->for_all_conversionless_edits(x.second, [&](VertexPair uv) {
                        m_costs_remaining[uv] -= current_min_cost;
                        return false;
                    });
                } else {
                    Q.emplace(current_min_cost, std::move(x.second));
                }
            }




            return lower_bound;
        }
    };
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_GREEDYWEIGHTEDPACKING_H
