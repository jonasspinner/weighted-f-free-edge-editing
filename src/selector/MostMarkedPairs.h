#ifndef WEIGHTED_F_FREE_EDGE_EDITING_MOSTMARKEDPAIRS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_MOSTMARKEDPAIRS_H


#include "../forbidden_subgraphs/SubgraphC4P4.h"


namespace selector {

    template<Options::FSG SetOfForbiddenSubgraphs>
    class MostMarkedPairs : public SelectorI {
    private:
        using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
        using Finder = typename Subgraph::Finder;

        const Graph &m_graph;
        const VertexPairMap<bool> &m_marked;
        const SubgraphStats<SetOfForbiddenSubgraphs> &m_subgraph_stats;

        Finder finder;
    public:
        MostMarkedPairs(const Graph &graph, const VertexPairMap<bool> &marked,
                        const SubgraphStats<SetOfForbiddenSubgraphs> &subgraph_stats) :
                m_graph(graph), m_marked(marked), m_subgraph_stats(subgraph_stats) {}

        Problem select_problem(Cost /*k*/) override {
            std::optional<Subgraph> max_subgraph;
            int max_num_marked_pairs = std::numeric_limits<int>::min();

            finder.find(m_graph, [&](Subgraph subgraph) {
                int num_marked_pairs = 0;
                for (VertexPair uv : subgraph.vertex_pairs())
                    if (m_marked[uv])
                        ++num_marked_pairs;

                if (num_marked_pairs > max_num_marked_pairs) {
                    max_num_marked_pairs = num_marked_pairs;
                    max_subgraph = std::move(subgraph);
                }
                return false;
            });

            if (!max_subgraph)
                return {{}, true};

            std::vector<VertexPair> pairs;
            for (VertexPair uv : max_subgraph->vertex_pairs())
                if (!m_marked[uv])
                    pairs.push_back(uv);

            std::sort(pairs.begin(), pairs.end(), [&](const VertexPair &uv, const VertexPair &xy) {
                return m_subgraph_stats.subgraphCount(uv) > m_subgraph_stats.subgraphCount(xy);
            });

            return {pairs, false};
        }
    };

    template
    class MostMarkedPairs<Options::FSG::C4P4>;
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_MOSTMARKEDPAIRS_H
