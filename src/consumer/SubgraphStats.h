//
// Created by jonas on 05.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHSTATS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHSTATS_H


#include "../interfaces/ConsumerI.h"
#include "../graph/VertexPairMap.h"
#include "../Instance.h"

class SubgraphStats : public ConsumerI {
private:
    VertexPairMap<size_t> subgraph_count_per_vertex_pair;
    size_t subgraph_count_per_vertex_pair_sum;
    size_t subgraph_count;
    std::vector<size_t> before_mark_subgraph_count;

    const Graph &graph;
    const VertexPairMap<bool> &m_forbidden;

    class State : public StateI {
        std::unique_ptr<StateI> copy() override {
            return std::make_unique<State>(*this);
        }
    };

public:
    SubgraphStats(std::shared_ptr<FinderI> finder_ptr, const Instance &instance, const VertexPairMap<bool> &forbidden)
            : ConsumerI(std::move(finder_ptr)), subgraph_count_per_vertex_pair(instance.graph.size()),
              subgraph_count_per_vertex_pair_sum(0), subgraph_count(0), graph(instance.graph), m_forbidden(forbidden) {}

    std::unique_ptr<StateI> initialize(Cost /*k*/) override {
        subgraph_count_per_vertex_pair = VertexPairMap<size_t>(graph.size());
        subgraph_count_per_vertex_pair_sum = 0;
        subgraph_count = 0;

        finder->find([&](Subgraph &&subgraph) {
            register_subgraph(subgraph);
            return false;
        });

        return std::make_unique<State>();
    }

    void before_mark_and_edit(StateI &/*state*/, VertexPair /*uv*/) override {}

    void after_mark_and_edit(StateI &/*state*/, VertexPair /*uv*/) override {}

    void before_mark(StateI &/*state*/, VertexPair /*uv*/) override {}

    void after_mark(StateI &/*state*/, VertexPair uv) override {
        subgraph_count_per_vertex_pair_sum -= subgraph_count_per_vertex_pair[uv];
        before_mark_subgraph_count.push_back(subgraph_count_per_vertex_pair[uv]);
        subgraph_count_per_vertex_pair[uv] = 0;

        verify();
    }

    void before_edit(StateI &/*state*/, VertexPair uv) override {
        assert(m_forbidden[uv]);
        finder->find_near(uv, [&](Subgraph &&subgraph) {
            remove_subgraph(subgraph);
            return false;
        });
        assert(subgraph_count_per_vertex_pair[uv] == 0);
    }

    void after_edit(StateI &/*state*/, VertexPair uv) override {
        finder->find_near(uv, [&](Subgraph &&subgraph) {
            register_subgraph(subgraph);
            return false;
        });
        assert(subgraph_count_per_vertex_pair[uv] == 0);
    }

    void after_unmark(StateI &/*state*/, VertexPair uv) override {
        subgraph_count_per_vertex_pair[uv] = before_mark_subgraph_count.back();
        subgraph_count_per_vertex_pair_sum += before_mark_subgraph_count.back();
        before_mark_subgraph_count.pop_back();
        verify();
    }

private:
    void register_subgraph(const Subgraph &subgraph) {
        subgraph_count++;
        for (VertexPair uv : subgraph.vertexPairs()) {
            if (!m_forbidden[uv]) {
                subgraph_count_per_vertex_pair[uv]++;
                subgraph_count_per_vertex_pair_sum++;
            }
        }
    }

    void remove_subgraph(const Subgraph &subgraph) {
        subgraph_count--;
        for (VertexPair uv : subgraph.vertexPairs()) {
            if (!m_forbidden[uv]) {
                subgraph_count_per_vertex_pair[uv]--;
                subgraph_count_per_vertex_pair_sum--;
            }
        }
    }

    void verify() {
#ifndef NDEBUG
        VertexPairMap<size_t> debug_sg_per_vertex_pair(graph.size());
        size_t debug_sg_count = 0;

        finder->find([&](Subgraph &&subgraph) {
            debug_sg_count++;

            for (VertexPair uv : subgraph.vertexPairs()) {
                if (!m_forbidden[uv]) {
                    debug_sg_per_vertex_pair[uv]++;
                }
            }
            return false;
        });

        assert(debug_sg_count == subgraph_count);

        for (VertexPair uv : graph.vertexPairs()) {
            assert(debug_sg_per_vertex_pair[uv] == subgraph_count_per_vertex_pair[uv]);
            assert(!m_forbidden[uv] || debug_sg_per_vertex_pair[uv]);
        }
#endif
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHSTATS_H
