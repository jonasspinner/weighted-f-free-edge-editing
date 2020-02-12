//
// Created by jonas on 12.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIRPREMARKER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIRPREMARKER_H


class VertexPairPreMarker {
    std::vector <VertexPair> m_vertex_pairs;
    std::vector <size_t> m_indices;
    std::vector<bool> m_needs_reset;

    const VertexPairMap <Cost> &m_costs;
    const std::vector<ConsumerI *> &m_consumers;

    VertexPairMap<bool> &m_marked;

public:
    VertexPairPreMarker(const Graph &graph, const VertexPairMap <Cost> &costs,
                        const std::vector<ConsumerI *> &consumers, VertexPairMap<bool> &marked) :
            m_indices({0}), m_costs(costs), m_consumers(consumers), m_marked(marked) {
        m_vertex_pairs.reserve(graph.size() * (graph.size() - 1) / 2);
        for (VertexPair uv : graph.vertexPairs())
            m_vertex_pairs.push_back(uv);

        std::sort(m_vertex_pairs.begin(), m_vertex_pairs.end(),
                  [&](VertexPair uv, VertexPair xy) { return m_costs[uv] > m_costs[xy]; });

        m_needs_reset.resize(m_vertex_pairs.size());
    }

    void mark(Cost k) {
        // All vertex pairs would be marked. Skip the work as only one selector call is necessary.
        if (k < m_costs[m_vertex_pairs.back()]) {
            m_indices.push_back(m_indices.back());
            return;
        };

        size_t index = m_indices.back();
        while (index < m_vertex_pairs.size() && m_costs[m_vertex_pairs[index]] > k) {
            VertexPair uv = m_vertex_pairs[index];
            if (!m_marked[uv]) {
                for (auto &c : m_consumers) c->before_mark(uv);
                m_marked[uv] = true;
                for (auto &c : m_consumers) c->after_mark(uv);
                m_needs_reset[index] = true;
            }
            ++index;
        }
        m_indices.push_back(index);
    }

    void unmark() {
        size_t end = m_indices.back();
        m_indices.pop_back();
        size_t begin = m_indices.back();
        for (size_t index = end - 1; begin <= index && index < end; --index) {
            if (m_needs_reset[index]) {
                auto uv = m_vertex_pairs[index];
                assert(m_marked[uv]);

                m_marked[uv] = false;
                for (auto &c : m_consumers) c->after_unmark(uv);

                m_needs_reset[index] = false;
            }
        }
    }

};

class PreMarkerGuard {
    VertexPairPreMarker &m_marker;
    bool m_marked = false;
public:
    explicit PreMarkerGuard(VertexPairPreMarker &marker) : m_marker(marker) {}

    void mark(Cost k) {
        assert(!m_marked);
        m_marker.mark(k);
        m_marked = true;
    }

    ~PreMarkerGuard() { if (m_marked) m_marker.unmark(); }
};

#endif //WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIRPREMARKER_H
