#ifndef WEIGHTED_F_FREE_EDGE_EDITING_EDITSTATE_H
#define WEIGHTED_F_FREE_EDGE_EDITING_EDITSTATE_H

#include "../graph/Graph.h"
#include "../graph/VertexPairMap.h"


class EditState {
    Graph m_graph;
    VertexPairMap<Cost> m_costs;
    VertexPairMap<bool> m_edited;
    VertexPairMap<bool> m_marked;

    std::vector<VertexPair> m_edits;

public:
    EditState(Graph graph, VertexPairMap<Cost> costs) :
            m_graph(std::move(graph)), m_costs(std::move(costs)), m_edited(m_graph.size()), m_marked(m_graph.size()) {}

    [[nodiscard]] constexpr const Graph &graph() const {
        return m_graph;
    }

    [[nodiscard]] constexpr const VertexPairMap<Cost> &cost_map() const {
        return m_costs;
    }

    [[nodiscard]] constexpr const VertexPairMap<bool> &marked_map() const {
        return m_marked;
    }

    [[nodiscard]] constexpr const std::vector<VertexPair> &edits() const {
        return m_edits;
    }

    [[nodiscard]] Cost cost(VertexPair uv) const {
        return m_costs[uv];
    }

    [[nodiscard]] bool is_edited(VertexPair uv) const {
        return m_edited[uv];
    }

    [[nodiscard]] bool is_marked(VertexPair uv) const {
        return m_marked[uv];
    }

    void edit_edge(VertexPair uv) {
        assert(is_marked(uv));
        assert(!is_edited(uv));

        m_graph.toggleEdge(uv);
        m_edits.push_back(uv);
        m_edited[uv] = true;
    }

    void unedit_edge(VertexPair uv) {
        assert(is_marked(uv));
        assert(is_edited(uv));

        m_graph.toggleEdge(uv);
        m_edits.pop_back();
        m_edited[uv] = false;
    }

    void mark_edge(VertexPair uv) {
        assert(!is_marked(uv));

        m_marked[uv] = true;
    }

    void unmark_edge(VertexPair uv) {
        assert(is_marked(uv));

        m_marked[uv] = false;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_EDITSTATE_H
