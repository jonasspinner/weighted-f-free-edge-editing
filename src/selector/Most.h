//
// Created by jonas on 18.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_MOST_H
#define WEIGHTED_F_FREE_EDGE_EDITING_MOST_H


#include "../interfaces/SelectorI.h"
#include "../consumer/SubgraphStats.h"


class Most : public SelectorI {
private:
    const VertexPairMap<bool> &m_marked;
    const SubgraphStats &m_subgraph_stats;
    Graph m_used;

public:
    Most(std::shared_ptr<FinderI> finder_ptr, const VertexPairMap<bool> &marked, const SubgraphStats &subgraph_stats) : SelectorI(
            std::move(finder_ptr)), m_marked(marked), m_subgraph_stats(subgraph_stats), m_used(m_marked.size()) {}

    Problem result(Cost /*k*/) override {

        // TODO: This selector triggers the assertion in the loop over the problem pairs: assert(!m_marked[uv]);

        size_t max_subgraph_count = 0;
        std::vector<VertexPair> pairs;
        for (VertexPair uv : Graph::VertexPairs(m_marked.size())) {
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
            finder->find_near(uv, m_used, [&](Subgraph &&subgraph) {

                for (VertexPair xy : subgraph.vertexPairs())
                    if (!m_marked[xy])
                        current_pairs.emplace_back(m_subgraph_stats.subgraphCount(xy), xy);

                std::sort(current_pairs.begin(), current_pairs.end(), std::greater<>());

                if (best_pairs.empty()) {
                    best_pairs = std::move(current_pairs);
                } else {
                    auto bend = best_pairs.end();
                    auto cend = current_pairs.end();

                    auto [bit, cit] = std::mismatch(best_pairs.begin(), bend, current_pairs.begin(), cend,
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

        for (auto [_, uv] : best_pairs) {
            assert(!m_marked[uv]);
            problem.pairs.push_back(uv);
        }

        return problem;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_MOST_H
