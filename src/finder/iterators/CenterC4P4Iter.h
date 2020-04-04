//
// Created by jonas on 11.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_CENTERC4P4ITER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_CENTERC4P4ITER_H

#include "../../graph/Subgraph.h"


class CenterC4P4Iter {
private:
    template <typename F, typename G, typename H, typename I>
    class Iterator {
    private:
        const Graph &m_graph;
        // Graph::Vertices::Iterator u_iter, v_iter, w_iter, x_iter;
        // Graph::AdjRow V, A, B;
        Graph::Vertices::Iterator u_iter;
        Graph::AdjRow V, A, B;

        Graph::RowVertices::Iterator v_iter = Graph::iterate(V).begin();
        Graph::RowVertices::Iterator a_iter = Graph::iterate(A).begin();
        Graph::RowVertices::Iterator b_iter = Graph::iterate(B).begin();

        Vertex u, v, a, b;

        Subgraph subgraph{};

        enum class Pos {
            start, return_1, return_2, end
        } pos;

        F neighbors; G non_neighbors; H valid_edge; I valid_non_edge;
    public:
        using value_type = Subgraph;
        using difference_type = std::ptrdiff_t;
        using pointer = const value_type *;
        using reference = const value_type &;
        using iterator_category = std::forward_iterator_tag;

        explicit Iterator(const Graph &graph, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge, Pos pos = Pos::start) :
            m_graph(graph),
            u_iter(m_graph.vertices().begin()),
            v_iter(Graph::iterate(V).begin()),
            a_iter(Graph::iterate(A).begin()),
            b_iter(Graph::iterate(B).begin()),
            u(0), v(0), a(0), b(0), pos(pos),
            neighbors(neighbors), non_neighbors(non_neighbors), valid_edge(valid_edge), valid_non_edge(valid_non_edge) {}

        Subgraph operator*() const {
            return subgraph;
        }

        Iterator &operator++() {
            switch (pos) {
                case Pos::start:
                    goto start_label;
                case Pos::return_1:
                    goto return_1_label;
                case Pos::return_2:
                    goto return_2_label;
                case Pos::end:
                    goto end_label;
            }

            start_label:
            for (; u_iter != m_graph.vertices().end(); ++u_iter) {
                u = *u_iter;
                V = neighbors(u);


                for (v_iter = Graph::iterate(V).begin(); v_iter != Graph::iterate(V).end(); ++v_iter) {
                    v = *v_iter;

                    A = neighbors(u) & non_neighbors(v);
                    B = neighbors(v) & non_neighbors(u);

                    for (a_iter = Graph::iterate(A).begin(); a_iter != Graph::iterate(A).end(); ++a_iter) {
                        a = *a_iter;

                        for (b_iter = Graph::iterate(B).begin(); b_iter != Graph::iterate(B).end();) {
                            b = *b_iter;

                            if (valid_non_edge({a, b}) && a < b) {
                                // P_4
                                pos = Pos::return_1;

                                subgraph = {a, u, v, b};
                                return *this;
                            }
                            return_1_label:
                            if (valid_edge({a, b}) && a < std::min({u, v, b}) && u < b) {
                                // C_4
                                subgraph = {a, u, v, b};
                                return *this;
                            }
                            return_2_label:
                            ++b_iter;
                        }
                    }
                }
            }
            pos = Pos::end;
            end_label:
            return *this;
        }

        bool operator==(const Iterator &other) const { return pos = Pos::end; }

        bool operator!=(const Iterator &other) const { return !(*this == other); }
    };

    template <typename F, typename G, typename H, typename I>
    class ForbiddenSubgraphs {
        const Graph &m_graph;
        F neighbors; G non_neighbors; H valid_edge; I valid_non_edge;
    public:
        ForbiddenSubgraphs(const Graph &graph, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) : m_graph(graph), neighbors(neighbors), non_neighbors(non_neighbors), valid_edge(valid_edge), valid_non_edge(valid_non_edge) {}
        [[nodiscard]] Iterator<F, G, H, I> begin() const { return Iterator(m_graph, neighbors, non_neighbors, valid_edge, valid_non_edge); }
        [[nodiscard]] Iterator<F, G, H, I> end() const { return Iterator(m_graph, neighbors, non_neighbors, valid_edge, valid_non_edge, Iterator<F, G, H, I>::Pos::end); }
    };
public:
    CenterC4P4Iter() = default;

    static inline auto neighbors(const Graph &graph) {
        return [&](Vertex u) -> const Graph::AdjRow& { return  graph.m_adj[u]; };
    }

    static inline auto neighbors(const Graph &graph, const Graph &forbidden) {
        return [&](Vertex u) { return graph.m_adj[u] - forbidden.m_adj[u]; };
    }

    static inline auto non_neighbors(const Graph &graph) {
        return [&](Vertex u) { auto result = ~graph.m_adj[u]; result[u] = false; return result; };
    }

    static inline auto non_neighbors(const Graph &graph, const Graph &forbidden) {
        return [&](Vertex u) { auto result = ~graph.m_adj[u] - forbidden.m_adj[u]; result[u] = false; return result; };
    }

    static inline auto valid_edge(const Graph &graph) {
        return [&](VertexPair xy) { return graph.hasEdge(xy); };
    }

    static inline auto valid_edge(const Graph &graph, const Graph &forbidden) {
        return [&](VertexPair xy) { return graph.hasEdge(xy) && !forbidden.hasEdge(xy); };
    }

    static inline auto valid_non_edge(const Graph &graph) {
        return [&](VertexPair xy) { return !graph.hasEdge(xy); };
    }

    static inline auto valid_non_edge(const Graph &graph, const Graph &forbidden) {
        return [&](VertexPair xy) { return !graph.hasEdge(xy) && !forbidden.hasEdge(xy); };
    }

    static auto forbidden_subgraphs(const Graph &graph) {
        ForbiddenSubgraphs(graph, neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
    }

    static auto forbidden_subgraphs(const Graph &graph, const Graph &forbidden) {
        ForbiddenSubgraphs(graph, neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

private:
//    template<typename F, typename G, typename H, typename I>
//    ForbiddenSubgraphs<F, G, H, I> forbidden_subgraphs(const Graph& graph, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
//        return ForbiddenSubgraphs(graph, neighbors, non_neighbors, valid_edge, valid_non_edge);
//    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_CENTERC4P4ITER_H
