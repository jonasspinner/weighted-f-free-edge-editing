#ifndef WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKING_H
#define WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKING_H


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

#endif //WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKING_H
