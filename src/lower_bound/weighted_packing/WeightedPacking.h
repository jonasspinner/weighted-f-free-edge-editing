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
     *
     * The weighted packing is maximal iff no forbidden subgraph can be inserted. This is exactly the case when all
     * forbidden subgraphs in the graph have at least one vertex pair which does not lead to a conversion and is
     * depleted. In the context of finders this is when find_with_duplicates finds no forbidden subgraph in m_graph
     * with the vertex pairs of m_depleted_graph marked as forbidden.
     */

    const Graph &m_graph;
    const VertexPairMap<Cost> &m_costs;
    const VertexPairMap<bool> &m_marked;
    const SubgraphStats &m_subgraph_stats;

    VertexPairMap<Cost> m_potential;
    Graph m_depleted_graph;
    Cost m_total_cost = 0;
    std::unique_ptr<FinderI> m_finder;

public:
    WeightedPacking(const Instance &instance, const VertexPairMap<bool> &marked, const SubgraphStats &subgraph_stats)
            : m_graph(instance.graph), m_costs(instance.costs), m_marked(marked), m_subgraph_stats(subgraph_stats),
              m_potential(m_graph.size()), m_depleted_graph(m_graph.size()) {}

    void clear() {
        for (VertexPair uv : m_graph.vertexPairs()) {
            m_potential[uv] = m_costs[uv];
        }
        m_depleted_graph.clearEdges();
        m_total_cost = 0;
    }

    /**
     * Modifies m_potential and m_depleted_graph.
     * @param subgraph
     * @param cost
     */
    void subtract_from_potential(const Subgraph &subgraph, Cost cost) {
        assert(cost > 0);
        assert(cost != invalid_cost);

        m_total_cost += cost;

        m_finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
            assert(!m_depleted_graph.hasEdge(uv));

            assert(cost <= m_potential[uv]);
            m_potential[uv] -= cost;
            if (m_potential[uv] == 0) {
                m_depleted_graph.setEdge(uv);
            }
            return false;
        });
    }

    /**
     * Modifies m_potential and m_depleted_graph.
     * @param subgraph
     * @param cost
     */
    void add_to_potential(const Subgraph &subgraph, Cost cost) {
        assert(cost > 0);
        assert(cost != invalid_cost);

        m_total_cost -= cost;

        m_finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
            m_potential[uv] += cost;
            m_depleted_graph.clearEdge(uv);
            return false;
        });
    }

    [[nodiscard]] Cost cost() const {
        return m_total_cost;
    }

    [[nodiscard]] bool is_depleted(VertexPair uv) const {
        return m_depleted_graph.hasEdge(uv);
    }

#ifndef NDEBUG

    bool is_maximal() const {
        return m_finder->find(m_graph, m_depleted_graph, [](auto) {
            return true;
        });
    }

#endif

    [[nodiscard]] bool is_valid() const {
        bool valid = true;
        for (VertexPair uv : Graph::VertexPairs(m_potential.size())) {
            valid |= 0 > m_potential[uv];
            assert(0 <= m_potential[uv]);
            valid |= m_potential[uv] > m_costs[uv];
            assert(m_potential[uv] <= m_costs[uv]);
        }

        for (VertexPair uv : Graph::VertexPairs(m_potential.size())) {
            valid |= m_potential[uv] == 0 != m_depleted_graph.hasEdge(uv);
            assert((m_potential[uv] == 0) == m_depleted_graph.hasEdge(uv));
        }
    }

    /**
     * Precondition:
     *      All pairs are not depleted.
     * Invariant:
     *      The depleted graph is the same after the call.
     *
     * @param finder
     * @param pairs
     * @param graph
     * @param depleted_graph
     * @param subgraph_stats
     * @return
     */
    std::pair<std::vector<Subgraph>, std::vector<size_t>>
    get_neighbors(const std::vector<VertexPair> &pairs) {

#ifndef NDEBUG
        for (auto uv : pairs) {
            assert(!m_depleted_graph.hasEdge(uv));
        }
#endif

        std::vector<Subgraph> candidates;
        std::vector<size_t> border(pairs.size() + 1);

        for (size_t i = 0; i < pairs.size(); ++i) {
            VertexPair uv = pairs[i];
            assert(!m_depleted_graph.hasEdge(uv));

            if (m_subgraph_stats.subgraphCount(uv) > 1) {
                m_finder->find_near_with_duplicates(uv, m_graph, m_depleted_graph, [&](Subgraph &&neighbor) {
#ifndef NDEBUG
                    m_finder->for_all_conversionless_edits(neighbor, [&](auto xy) {
                        assert(!m_depleted_graph.hasEdge(xy));
                        return false;
                    });
#endif
                    candidates.push_back(std::move(neighbor));
                    return false;
                });
            }
            border[i + 1] = candidates.size();

            // prevent subgraphs including uv to be counted twice
            m_depleted_graph.setEdge(uv);
        }

        // reset bound_graph
        for (VertexPair uv : pairs) {
            assert(m_depleted_graph.hasEdge(uv));
            m_depleted_graph.clearEdge(uv);
        }

        return {std::move(candidates), std::move(border)};
    }

    std::vector<VertexPair> get_neighbor_pairs(const Subgraph &subgraph) {
        std::vector<VertexPair> pairs;
        m_finder->for_all_conversionless_edits(subgraph, [&](auto uv) {
            assert(!m_depleted_graph.hasEdge(uv));
            if (!m_marked[uv])
                pairs.push_back(uv);
            return false;
        });
        return pairs;
    }

    [[nodiscard]] Cost calculate_min_cost(const Subgraph &subgraph) const {
        return m_finder->calculate_min_cost(subgraph, m_marked, m_potential);
    }
};

#endif //WEIGHTED_F_FREE_EDGE_EDITING_WEIGHTEDPACKING_H
