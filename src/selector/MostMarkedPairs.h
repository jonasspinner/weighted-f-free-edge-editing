#ifndef WEIGHTED_F_FREE_EDGE_EDITING_MOSTMARKEDPAIRS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_MOSTMARKEDPAIRS_H


namespace Selector {
    class MostMarkedPairs : public SelectorI {
    private:
        const Graph &m_graph;
        const VertexPairMap<bool> &m_marked;
        const SubgraphStats &m_subgraph_stats;

        std::shared_ptr<FinderI> finder;
    public:
        MostMarkedPairs(std::shared_ptr <FinderI> finder_ptr, const Graph &graph, const VertexPairMap<bool> &marked,
            const SubgraphStats &subgraph_stats) :
                m_graph(graph), m_marked(marked), m_subgraph_stats(subgraph_stats), finder(std::move(finder_ptr)) {}

        Problem select_problem(Cost /*k*/) override {
            Subgraph max_subgraph{};
            int max_num_marked_pairs = std::numeric_limits<int>::min();


            bool solved = true;

            finder->find(m_graph, [&](Subgraph &&subgraph) {
                solved = false;

                int num_marked_pairs = 0;
                for (VertexPair uv : subgraph.vertexPairs())
                    if (m_marked[uv])
                        ++num_marked_pairs;

                if (num_marked_pairs > max_num_marked_pairs) {
                    max_num_marked_pairs = num_marked_pairs;
                    max_subgraph = std::move(subgraph);
                }
                return false;
            });


            std::vector<VertexPair> pairs;
            for (VertexPair uv : max_subgraph.vertexPairs())
                if (!m_marked[uv])
                    pairs.push_back(uv);

            std::sort(pairs.begin(), pairs.end(), [&](const VertexPair &uv, const VertexPair &xy) {
                return m_subgraph_stats.subgraphCount(uv) > m_subgraph_stats.subgraphCount(xy);
            });

            return {pairs, solved};
        }
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_MOSTMARKEDPAIRS_H
