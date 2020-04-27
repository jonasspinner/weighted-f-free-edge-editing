#ifndef WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKING_H
#define WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKING_H


class WeightedPacking {
    /**
     * This class manages a weighted packing. A forbidden subgraph which is in the packing, has a cost associated with
     * it.
     *
     *   B := \{ (S, c) | S \subseteq V, c > 0 \}
     *
     * The lower bound is \sum_{(S, c) \in B} c. For a valid lower bound, the some restrictions must be enforced. We use
     * \phi(uv) to denote a potential which symbolizes how much more cost can be carried by the vertex pair.
     *
     *   \forall uv:
     *      0 <= \phi(uv) <= cost(uv)
     *      \phi(uv) = cost(uv) - \sum_{(S, c) in B, uv \in edits(S)} c
     *
     * As a consequence, an empty packing B = \empty results in \phi(uv) = cost(uv) \forall uv.
     *
     * To prevent the listing of forbidden subgraphs, that cannot be inserted into the packing, vertex pairs with
     * \phi(uv) = 0 are marked as depleted and not listed.
     *
     *   \forall uv: \phi(uv) = 0  <=>  uv \in E_{depleted}
     *
     * The goal is to use local search techniques to find ways to improve the packing. The FPT-Algorithm marks and edits
     * vertex pairs locally. The weighted packing can be updated respectively to remain valid.
     */

    VertexPairMap <Cost> m_potential;
    Graph m_depleted_graph;
    Cost m_total_cost = 0;
    std::unique_ptr <FinderI> m_finder;
    const Graph &m_graph;
    const VertexPairMap<Cost> &m_costs;
    std::vector <std::pair<Subgraph, Cost>> m_subgraphs_in_packing;

    /**
     * Modifies m_potenital and m_depleted_graph.
     * @param subgraph
     * @param cost
     */
    void subtract_from_potential(const Subgraph &subgraph, Cost cost) {
        assert(cost > 0);
        assert(cost != invalid_cost);

        m_finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
            assert(!m_saturated_vertex_pairs.hasEdge(uv));

            assert(cost <= m_potential[uv])
            m_potential[uv] -= cost;
            if (m_potential[uv] == 0) {
                m_depleted_graph.setEdge(uv);
            }
            return false;
        });
    }

    /**
     * Modifies m_potenital and m_depleted_graph.
     * @param subgraph
     * @param cost
     */
    void add_to_potential(const Subgraph &subgraph, Cost cost) {
        assert(cost > 0);
        assert(cost != invalid_cost);

        m_finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
            m_potential[uv] += cost;
            m_saturated_vertex_pairs.clearEdge(uv);
            return false;
        });
    }

    void insert_into_packing(const Subgraph &subgraph) {
        auto cost = finder->
        subtract_from_potential(subgraph)
        m_total_cost += cost;

        m_subgraphs_in_packing.emplace_back(subgraph, cost);
    }

    void remove_from_packing(size_t index) {
        const auto&[subgraph, cost] = m_subgraphs_in_packing[index];

        m_finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
            m_potential[uv] += cost;
            if (m_potential[uv] > 0) {
                m_saturated_vertex_pairs.clearEdge(uv);
            }
            return false;
        });
    }

    std::vector <Subgraph> neighbors(const Subgraph &subgraph) {
        std::vector <Subgraph> result;

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

        max = (c_S, S, _);

        for (size_t a_index = 0; a_index < n.size(); ++a_index) {
            c_a = insert_into_packing(a);
            for (size_t b_index = a_index + 1; b_index < n.size(); ++b_index) {
                c_b = insert_into_packing(b);
                if (both fit into packing)
                    max = (c_a + c_b, a, b);
                remove_from_packing(b, c_b);
            }
            remove_from_packing(a, c_a);
        }

        (c, a, b) = max;
        insert_into_packing(a);
        insert_into_packing(b);

    }

#ifndef NDEBUG
    bool is_maximal() {
        return m_finder->find(m_graph, m_saturated_vertex_pairs, [](auto) {
            return true;
        });
    }
#endif

    void local_search() {
        size_t max_iter = 5;
        size_t num_rounds_no_improvement = 0;

        for (size_t iter = 0; iter < max_iter; ++iter) {
            bool found_improvement = false;
            for (size_t index = 0; index < m_subgraphs_in_packing.size(); ++index) {
                found_improvement |= find_one_two_improvement(index);
            }
            if (!found_improvement)
                num_rounds_no_improvement++;
        }
    }

    void greedy_initialize() {
        ...
        assert(is_maximal());
    }

    Cost cost() const {
        return m_total_cost;
    }

    bool is_valid() const {
        bool valid = true;
        for (VertexPair uv : Graph::VertexPairs(m_potential.size())) {
            valid |= !(0 <= m_potential[uv]); assert(0 <= m_potential[uv]);
            valid |= !(m_potential[uv] <= m_costs[uv]); assert(m_potential[uv] <= m_costs[uv]);
        }

        auto other_potential = m_costs;
        for (const auto&[subgraph, cost] : m_subgraphs_in_packing) {
            m_finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
                other_potential[uv] -= cost;
                return false;
            });
        }

        for (VertexPair uv : Graph::VertexPairs(m_potential.size())) {
            valid |= !(other_potential[uv] == m_potential[uv]); assert(other_potential[uv] == m_potential[uv]);
        }

        for (VertexPair uv : Graph::VertexPairs(m_potential.size())) {
            valid |= !((m_potential[uv] == 0) == m_depleted_graph.hasEdge(uv)); assert((m_potential[uv] == 0) == m_depleted_graph.hasEdge(uv));
        }
    }

    void local_fix_marked(VertexPair uv) {

    }

    void local_fix_edited(VertexPair uv) {

    }
};

#endif //WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKING_H
