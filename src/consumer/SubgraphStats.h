#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHSTATS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHSTATS_H


#include "../consumer/ConsumerI.h"
#include "../graph/VertexPairMap.h"
#include "../Instance.h"
#include "../forbidden_subgraphs/SubgraphC4P4.h"
#include "../editor/EditState.h"


template <Options::FSG SetOfForbiddenSubgraphs>
class SubgraphStats final : public ConsumerI {
private:
    using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
    using Finder = typename Subgraph::Finder;

    Finder m_finder;

    const EditState *m_edit_state;

    VertexPairMap<size_t> m_subgraph_count_per_vertex_pair;
    size_t m_subgraph_count_per_vertex_pair_sum;
    size_t m_subgraph_count;
    std::vector<size_t> m_before_mark_subgraph_count;

    // TODO: Remove dependency on empty graph, i.e. build find_near_unique variant without forbidden graph.
    Graph m_empty_graph;

public:
    SubgraphStats(const EditState *edit_state)
            : m_edit_state(std::move(edit_state)), m_subgraph_count_per_vertex_pair(m_edit_state->graph().size()),
              m_subgraph_count_per_vertex_pair_sum(0), m_subgraph_count(0),
              m_empty_graph(m_edit_state->graph().size()) {}

    [[nodiscard]] constexpr size_t subgraph_count() const {
        return m_subgraph_count;
    }

    [[nodiscard]] constexpr size_t subgraph_count(VertexPair uv) const {
        return m_subgraph_count_per_vertex_pair[uv];
    }

    void initialize(Cost /*k*/) override {
        m_subgraph_count_per_vertex_pair = VertexPairMap<size_t>(m_edit_state->graph().size());
        m_subgraph_count_per_vertex_pair_sum = 0;
        m_subgraph_count = 0;

        m_finder.find(m_edit_state->graph(), [&](Subgraph subgraph) {
            register_subgraph(subgraph);
            return subgraph_iterators::IterationControl::Continue;
        });

        verify();
    }

    void remove_near_subgraphs(VertexPair uv) {
        assert(m_edit_state->is_marked(uv));
        verify();
        m_finder.find_near(uv, m_edit_state->graph(), m_empty_graph, [&](Subgraph subgraph) {
            remove_subgraph(subgraph);
            return subgraph_iterators::IterationControl::Continue;
        });
        assert(m_subgraph_count_per_vertex_pair[uv] == 0);
    }

    void register_near_subgraphs(VertexPair uv) {
        m_finder.find_near(uv, m_edit_state->graph(), m_empty_graph, [&](Subgraph subgraph) {
            register_subgraph(subgraph);
            return subgraph_iterators::IterationControl::Continue;
        });
        verify();
        assert(m_subgraph_count_per_vertex_pair[uv] == 0);
    }

    void before_edit(VertexPair uv) override {
        remove_near_subgraphs(uv);
    }

    void after_edit(VertexPair uv) override {
        register_near_subgraphs(uv);
    }

    void before_unedit(VertexPair uv) override {
        remove_near_subgraphs(uv);
    }

    void after_unedit(VertexPair uv) override {
        register_near_subgraphs(uv);
    }

    void after_mark(VertexPair uv) override {
        m_subgraph_count_per_vertex_pair_sum -= m_subgraph_count_per_vertex_pair[uv];
        m_before_mark_subgraph_count.push_back(m_subgraph_count_per_vertex_pair[uv]);
        m_subgraph_count_per_vertex_pair[uv] = 0;
        verify();
    }

    void after_unmark(VertexPair uv) override {
        m_subgraph_count_per_vertex_pair[uv] = m_before_mark_subgraph_count.back();
        m_subgraph_count_per_vertex_pair_sum += m_before_mark_subgraph_count.back();
        m_before_mark_subgraph_count.pop_back();
        verify();
    }

private:
    void register_subgraph(const Subgraph &subgraph) {
        m_subgraph_count++;
        for (VertexPair uv : subgraph.non_converting_edits()) {
            if (!m_edit_state->is_marked(uv)) {
                m_subgraph_count_per_vertex_pair[uv]++;
                m_subgraph_count_per_vertex_pair_sum++;
            }
        }
    }

    void remove_subgraph(const Subgraph &subgraph) {
        m_subgraph_count--;
        for (VertexPair uv : subgraph.non_converting_edits()) {
            if (!m_edit_state->is_marked(uv)) {
                m_subgraph_count_per_vertex_pair[uv]--;
                m_subgraph_count_per_vertex_pair_sum--;
            }
        }
    }

    void verify() {
#ifndef NDEBUG
        VertexPairMap<size_t> debug_sg_per_vertex_pair(m_edit_state->graph().size());
        size_t debug_sg_count = 0;

        m_finder.find(m_edit_state->graph(), [&](Subgraph subgraph) {
            debug_sg_count++;

            for (VertexPair uv : subgraph.non_converting_edits()) {
                if (!m_edit_state->is_marked(uv)) {
                    debug_sg_per_vertex_pair[uv]++;
                }
            }
            return subgraph_iterators::IterationControl::Continue;
        });

        assert(debug_sg_count == m_subgraph_count);

        for (VertexPair uv : m_edit_state->graph().vertex_pairs()) {
            assert(debug_sg_per_vertex_pair[uv] == m_subgraph_count_per_vertex_pair[uv]);
            assert(!m_edit_state->is_marked(uv) || debug_sg_per_vertex_pair[uv] == 0);
        }
#endif
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHSTATS_H
