#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHSTATS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHSTATS_H


#include "../consumer/ConsumerI.h"
#include "../graph/VertexPairMap.h"
#include "../Instance.h"
#include "../forbidden_subgraphs/SubgraphC4P4.h"


template <Options::FSG SetOfForbiddenSubgraphs>
class SubgraphStats final : public ConsumerI {
private:
    using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
    using Finder = typename Subgraph::Finder;

    Finder m_finder;

    VertexPairMap<size_t> subgraph_count_per_vertex_pair;
    size_t subgraph_count_per_vertex_pair_sum;
    size_t subgraph_count;
    std::vector<size_t> before_mark_subgraph_count;

    const Graph &m_graph;
    const VertexPairMap<bool> &m_marked;

    // TODO: Remove dependency on empty graph, i.e. build find_near_unique variant without forbidden graph.
    Graph m_empty_graph;

public:
    SubgraphStats(const Instance &instance, const VertexPairMap<bool> &marked)
            : subgraph_count_per_vertex_pair(instance.graph.size()),
              subgraph_count_per_vertex_pair_sum(0), subgraph_count(0), m_graph(instance.graph), m_marked(marked),
              m_empty_graph(m_graph.size()) {}

    [[nodiscard]] constexpr size_t subgraphCount() const {
        return subgraph_count;
    }

    [[nodiscard]] constexpr size_t subgraphCount(VertexPair uv) const {
        return subgraph_count_per_vertex_pair[uv];
    }

    void initialize(Cost /*k*/) override {
        subgraph_count_per_vertex_pair = VertexPairMap<size_t>(m_graph.size());
        subgraph_count_per_vertex_pair_sum = 0;
        subgraph_count = 0;

        m_finder.find(m_graph, [&](Subgraph subgraph) {
            register_subgraph(subgraph);
            return subgraph_iterators::IterationControl::Continue;
        });

        verify();
    }

    void remove_near_subgraphs(VertexPair uv) {
        assert(m_marked[uv]);
        verify();
        m_finder.find_near(uv, m_graph, m_empty_graph, [&](Subgraph subgraph) {
            remove_subgraph(subgraph);
            return subgraph_iterators::IterationControl::Continue;
        });
        assert(subgraph_count_per_vertex_pair[uv] == 0);
    }

    void register_near_subgraphs(VertexPair uv) {
        m_finder.find_near(uv, m_graph, m_empty_graph, [&](Subgraph subgraph) {
            register_subgraph(subgraph);
            return subgraph_iterators::IterationControl::Continue;
        });
        verify();
        assert(subgraph_count_per_vertex_pair[uv] == 0);
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
        subgraph_count_per_vertex_pair_sum -= subgraph_count_per_vertex_pair[uv];
        before_mark_subgraph_count.push_back(subgraph_count_per_vertex_pair[uv]);
        subgraph_count_per_vertex_pair[uv] = 0;
        verify();
    }

    void after_unmark(VertexPair uv) override {
        subgraph_count_per_vertex_pair[uv] = before_mark_subgraph_count.back();
        subgraph_count_per_vertex_pair_sum += before_mark_subgraph_count.back();
        before_mark_subgraph_count.pop_back();
        verify();
    }

private:
    void register_subgraph(const Subgraph &subgraph) {
        subgraph_count++;
        for (VertexPair uv : subgraph.non_converting_edits()) {
            if (!m_marked[uv]) {
                subgraph_count_per_vertex_pair[uv]++;
                subgraph_count_per_vertex_pair_sum++;
            }
        }
    }

    void remove_subgraph(const Subgraph &subgraph) {
        subgraph_count--;
        for (VertexPair uv : subgraph.non_converting_edits()) {
            if (!m_marked[uv]) {
                subgraph_count_per_vertex_pair[uv]--;
                subgraph_count_per_vertex_pair_sum--;
            }
        }
    }

    void verify() {
#ifndef NDEBUG
        VertexPairMap<size_t> debug_sg_per_vertex_pair(m_graph.size());
        size_t debug_sg_count = 0;

        m_finder.find(m_graph, [&](Subgraph subgraph) {
            debug_sg_count++;

            for (VertexPair uv : subgraph.non_converting_edits()) {
                if (!m_marked[uv]) {
                    debug_sg_per_vertex_pair[uv]++;
                }
            }
            return subgraph_iterators::IterationControl::Continue;
        });

        assert(debug_sg_count == subgraph_count);

        for (VertexPair uv : m_graph.vertexPairs()) {
            assert(debug_sg_per_vertex_pair[uv] == subgraph_count_per_vertex_pair[uv]);
            assert(!m_marked[uv] || debug_sg_per_vertex_pair[uv] == 0);
        }
#endif
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHSTATS_H
