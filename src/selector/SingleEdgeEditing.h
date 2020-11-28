#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SINGLEEDGEEDITING_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SINGLEEDGEEDITING_H


#include "../forbidden_subgraphs/SubgraphC4P4.h"


namespace selector {

    template<Options::FSG SetOfForbiddenSubgraphs>
    class SingleEdgeEditing : public SelectorI {
    private:
        using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
        using Finder = typename Subgraph::Finder;

        const EditState *m_edit_state;
        const SubgraphStats<SetOfForbiddenSubgraphs> *m_subgraph_stats;

        Finder finder;
    public:
        SingleEdgeEditing(const EditState *edit_state,
                          const SubgraphStats<SetOfForbiddenSubgraphs> *subgraph_stats) :
                m_edit_state(edit_state), m_subgraph_stats(subgraph_stats) {}

        [[nodiscard]] RecursionType recursion_type() const override { return RecursionType::VertexPair; }

        Problem select_problem(Cost /*k*/) override {

            if (m_subgraph_stats->subgraph_count() == 0)
                return {{}, true};

            VertexPair max_pair{0, 1};
            size_t max_count = 0;

            for (VertexPair uv : m_edit_state->graph().vertex_pairs()) {
                auto count = m_subgraph_stats->subgraph_count(uv);
                if (!m_edit_state->is_marked(uv) && count > max_count) {
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

    template
    class SingleEdgeEditing<Options::FSG::C4P4>;
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_SINGLEEDGEEDITING_H
