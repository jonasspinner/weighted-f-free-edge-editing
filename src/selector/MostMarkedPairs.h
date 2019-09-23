//
// Created by jonas on 29.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_MOSTMARKEDPAIRS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_MOSTMARKEDPAIRS_H


namespace Selector {
    class MostMarkedPairs : public SelectorI {
    private:
        const VertexPairMap<bool> &m_marked;
        const SubgraphStats &m_subgraph_stats;
    public:
        MostMarkedPairs(std::shared_ptr <FinderI> finder_ptr, const VertexPairMap<bool> &marked, const SubgraphStats &subgraph_stats) : SelectorI(
                std::move(finder_ptr)), m_marked(marked), m_subgraph_stats(subgraph_stats) {}

        Problem select_problem(Cost /*k*/) override {
            Subgraph min_subgraph{};
            size_t min_num_marked_pairs = std::numeric_limits<size_t>::max();


            bool solved = true;

            finder->find([&](Subgraph &&subgraph) {
                solved = false;

                size_t num_marked_pairs = 0;
                for (VertexPair uv : subgraph.vertexPairs())
                    if (m_marked[uv])
                        ++num_marked_pairs;

                if (num_marked_pairs < min_num_marked_pairs) {
                    min_num_marked_pairs = num_marked_pairs;
                    min_subgraph = std::move(subgraph);
                }
                return false;
            });


            std::vector<VertexPair> pairs;
            for (VertexPair uv : min_subgraph.vertexPairs())
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
