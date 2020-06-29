#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHSTATS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHSTATS_H


#include "../consumer/ConsumerI.h"
#include "../graph/VertexPairMap.h"
#include "../Instance.h"
#include "../forbidden_subgraphs/SubgraphC4P4.h"


class SubgraphStats : public ConsumerI {
public:
    [[nodiscard]] virtual size_t subgraphCount() const = 0;

    [[nodiscard]] virtual size_t subgraphCount(VertexPair uv) const = 0;
};


namespace subgraph_stats {

template <Options::FSG SetOfForbiddenSubgraphs>
class SubgraphStatsT : public SubgraphStats {
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

    Graph m_empty_graph;

public:
    SubgraphStatsT(const Instance &instance, const VertexPairMap<bool> &marked)
            : subgraph_count_per_vertex_pair(instance.graph.size()),
              subgraph_count_per_vertex_pair_sum(0), subgraph_count(0), m_graph(instance.graph), m_marked(marked),
              m_empty_graph(m_graph.size()) {}

    [[nodiscard]] size_t subgraphCount() const override {
        return subgraph_count;
    }

    [[nodiscard]] size_t subgraphCount(VertexPair uv) const override {
        return subgraph_count_per_vertex_pair[uv];
    }

    void initialize(Cost /*k*/) override {
        subgraph_count_per_vertex_pair = VertexPairMap<size_t>(m_graph.size());
        subgraph_count_per_vertex_pair_sum = 0;
        subgraph_count = 0;

        m_finder.find_unique(m_graph, [&](Subgraph subgraph) {
            register_subgraph(subgraph);
            return false;
        });

        verify();
    }

    void remove_near_subgraphs(VertexPair uv) {
        assert(m_marked[uv]);
        verify();
        m_finder.find_near_unique(uv, m_graph, m_empty_graph, [&](Subgraph subgraph) {
            remove_subgraph(subgraph);
            return false;
        });
        assert(subgraph_count_per_vertex_pair[uv] == 0);
    }

    void register_near_subgraphs(VertexPair uv) {
        m_finder.find_near_unique(uv, m_graph, m_empty_graph, [&](Subgraph subgraph) {
            register_subgraph(subgraph);
            return false;
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
        for (VertexPair uv : subgraph.vertex_pairs()) {
            if (!m_marked[uv]) {
                subgraph_count_per_vertex_pair[uv]++;
                subgraph_count_per_vertex_pair_sum++;
            }
        }
    }

    void remove_subgraph(const Subgraph &subgraph) {
        subgraph_count--;
        for (VertexPair uv : subgraph.vertex_pairs()) {
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

        m_finder.find_unique(m_graph, [&](Subgraph subgraph) {
            debug_sg_count++;

            for (VertexPair uv : subgraph.vertex_pairs()) {
                if (!m_marked[uv]) {
                    debug_sg_per_vertex_pair[uv]++;
                }
            }
            return false;
        });

        assert(debug_sg_count == subgraph_count);

        for (VertexPair uv : m_graph.vertexPairs()) {
            assert(debug_sg_per_vertex_pair[uv] == subgraph_count_per_vertex_pair[uv]);
            assert(!m_marked[uv] || debug_sg_per_vertex_pair[uv] == 0);
        }
#endif
    }
};


inline std::unique_ptr<SubgraphStats> make(Options::FSG fsg, const Instance &instance, const VertexPairMap<bool> &marked) {
    switch (fsg) {
        case Options::FSG::C4P4:
            return std::make_unique<SubgraphStatsT<Options::FSG::C4P4>>(instance, marked);
        default:
            throw std::runtime_error("SubgraphStatsT not specialized for given forbidden subgraphs.");
    }
}

}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHSTATS_H
