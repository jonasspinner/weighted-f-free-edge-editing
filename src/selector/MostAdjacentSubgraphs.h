#ifndef WEIGHTED_F_FREE_EDGE_EDITING_MOSTADJACENTSUBGRAPHS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_MOSTADJACENTSUBGRAPHS_H


#include "SelectorI.h"
#include "../consumer/SubgraphStats.h"
#include "../forbidden_subgraphs/SubgraphC4P4.h"


namespace selector {

    template<Options::FSG SetOfForbiddenSubgraphs>
    class MostAdjacentSubgraphs : public SelectorI {
    private:
        using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
        using Finder = typename Subgraph::Finder;

        const Graph &m_graph;
        const VertexPairMap<bool> &m_marked;
        const SubgraphStats<SetOfForbiddenSubgraphs> &m_subgraph_stats;
        Graph m_used;

        Finder m_finder;
    public:
        MostAdjacentSubgraphs(const Graph &graph, const VertexPairMap<bool> &marked,
                              const SubgraphStats<SetOfForbiddenSubgraphs> &subgraph_stats) :
                m_graph(graph), m_marked(marked), m_subgraph_stats(subgraph_stats),
                m_used(m_marked.size()) {}

        Problem select_problem(Cost /*k*/) override {

            size_t max_subgraph_count = 0;
            std::vector<VertexPair> pairs;
            for (VertexPair uv : m_graph.vertexPairs()) {
                size_t subgraph_count = m_subgraph_stats.subgraphCount(uv);
                if (subgraph_count > max_subgraph_count) {
                    max_subgraph_count = subgraph_count;
                    pairs = {uv};
                } else if (subgraph_count == max_subgraph_count && (subgraph_count > 1 || pairs.empty())) {
                    pairs.push_back(uv);
                }
            }

            std::vector<std::pair<size_t, VertexPair>> best_pairs, current_pairs;
            for (VertexPair uv : pairs) {
                m_finder.find_near(uv, m_graph, m_used, [&](const Subgraph &subgraph) {
                    current_pairs.clear();

                    for (VertexPair xy : subgraph.vertex_pairs())
                        if (!m_marked[xy])
                            current_pairs.emplace_back(m_subgraph_stats.subgraphCount(xy), xy);

                    std::sort(current_pairs.begin(), current_pairs.end(), std::greater<>());

                    if (best_pairs.empty()) {
                        best_pairs = std::move(current_pairs);
                    } else {
                        auto bend = best_pairs.end();
                        auto cend = current_pairs.end();

                        auto[bit, cit] = std::mismatch(best_pairs.begin(), bend, current_pairs.begin(), cend,
                                                       [](const auto &a, const auto &b) { return a.first == b.first; });

                        if (cit == cend || (bit != bend && *bit < *cit))
                            best_pairs = std::move(current_pairs);
                    }

                    return false;
                });

                m_used.setEdge(uv);
            }

            for (VertexPair uv : pairs)
                m_used.clearEdge(uv);


            Problem problem;
            problem.solved = (m_subgraph_stats.subgraphCount() == 0);

#ifndef NDEBUG
            if (problem.solved) {
                assert(!finder.find(m_graph, [](Subgraph) { return true; }));
            }
#endif

            for (auto[_, uv] : best_pairs) {
                assert(!m_marked[uv]);
                problem.pairs.push_back(uv);
            }

            return problem;
        }
    };

    template
    class MostAdjacentSubgraphs<Options::FSG::C4P4>;
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_MOSTADJACENTSUBGRAPHS_H
