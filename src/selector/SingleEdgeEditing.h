#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SINGLEEDGEEDITING_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SINGLEEDGEEDITING_H


namespace Selector {
    class SingleEdgeEditing : public SelectorI {
    private:
        const Graph &m_graph;
        const VertexPairMap<bool> &m_marked;
        const SubgraphStats &m_subgraph_stats;

        std::shared_ptr<FinderI> finder;
    public:
        SingleEdgeEditing(std::shared_ptr<FinderI> finder_ptr, const Graph &graph, const VertexPairMap<bool> &marked,
                          const SubgraphStats &subgraph_stats) :
                m_graph(graph), m_marked(marked), m_subgraph_stats(subgraph_stats),
                finder(std::move(finder_ptr)) {}

        [[nodiscard]] RecursionType recursion_type() const override { return RecursionType::VertexPair; }

        Problem select_problem(Cost /*k*/) override {

            if (m_subgraph_stats.subgraphCount() == 0)
                return {{}, true};

            VertexPair max_pair{0, 1};
            size_t max_count = 0;

            for (VertexPair uv : m_graph.vertexPairs()) {
                auto count = m_subgraph_stats.subgraphCount(uv);
                if (!m_marked[uv] && count > max_count) {
                    max_count = count;
                    max_pair = uv;
                }
            }

            if (max_count > 0) {
                return {{max_pair}, false};
            } else {
                // Every subgraph is fully marked.
                return {{}, false};
            }
        }
    };
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_SINGLEEDGEEDITING_H
