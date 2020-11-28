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

        const EditState *m_edit_state;
        const SubgraphStats<SetOfForbiddenSubgraphs> *m_subgraph_stats;
        Graph m_used;

        Finder m_finder;
    public:
        MostAdjacentSubgraphs(const EditState *edit_state,
                              const SubgraphStats<SetOfForbiddenSubgraphs> *subgraph_stats) :
                m_edit_state(edit_state), m_subgraph_stats(subgraph_stats),
                m_used(m_edit_state->graph().number_of_vertices()) {}

        Problem select_problem(Cost /*k*/) override {
            const auto &graph = m_edit_state->graph();

            size_t max_subgraph_count = 0;
            std::vector<VertexPair> pairs;
            for (VertexPair uv : graph.vertex_pairs()) {
                size_t subgraph_count = m_subgraph_stats->subgraph_count(uv);
                if (subgraph_count > max_subgraph_count) {
                    max_subgraph_count = subgraph_count;
                    pairs = {uv};
                } else if (subgraph_count == max_subgraph_count && (subgraph_count > 1 || pairs.empty())) {
                    pairs.push_back(uv);
                }
            }

            std::vector<std::pair<size_t, VertexPair>> best_pairs, current_pairs;
            for (VertexPair uv : pairs) {
                m_finder.find_near(uv, graph, m_used, [&](const Subgraph &subgraph) {
                    current_pairs.clear();

                    for (VertexPair xy : subgraph.vertex_pairs())
                        if (!m_edit_state->is_marked(xy))
                            current_pairs.emplace_back(m_subgraph_stats->subgraph_count(xy), xy);

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

                    return subgraph_iterators::IterationControl::Continue;
                });

                m_used.set_edge(uv);
            }

            for (VertexPair uv : pairs)
                m_used.reset_edge(uv);


            Problem problem;
            problem.solved = (m_subgraph_stats->subgraph_count() == 0);

#ifndef NDEBUG
            if (problem.solved) {
                assert(subgraph_iterators::IterationExit::Normal == m_finder.find(graph, [](Subgraph) { return subgraph_iterators::IterationControl::Break; }));
            }
#endif

            for (auto[_, uv] : best_pairs) {
                assert(!m_edit_state->is_marked(uv));
                problem.pairs.push_back(uv);
            }

            return problem;
        }
    };

    template
    class MostAdjacentSubgraphs<Options::FSG::C4P4>;
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_MOSTADJACENTSUBGRAPHS_H
