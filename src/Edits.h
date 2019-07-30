//
// Created by jonas on 26.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_EDITS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_EDITS_H


#include <vector>
#include "graph/Graph.h"

class Edits {
    std::vector<VertexPair> m_pairs;
    Graph m_graph;

public:
    explicit Edits(const Graph &graph) : m_graph(graph.n_vertices()) {}

    [[nodiscard]] bool is_marked(VertexPair uv) const {
        return m_graph.has_edge(uv);
    }

    [[nodiscard]] inline size_t size() const {
        return m_pairs.size();
    }

    [[nodiscard]] inline bool empty() const {
        return m_pairs.empty();
    }

    void edit_and_mark(VertexPair uv) {
        assert(!is_marked(uv));
        m_pairs.push_back(uv);
        m_graph.set_edge(uv);
    }

    [[nodiscard]] inline auto last_edit() const {
        return m_pairs.back();
    }

    void undo_last_edit() {
        assert(!empty());
        m_pairs.pop_back();
    }

    void unmark(VertexPair uv) {
        assert(is_marked(uv));
        m_graph.clear_edge(uv);
    }

    [[nodiscard]] const Graph &marked() const {
        return m_graph;
    }

    [[nodiscard]] const std::vector<VertexPair> &pairs() const {
        return m_pairs;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_EDITS_H
