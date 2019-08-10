//
// Created by jonas on 03.07.19.
//

#ifndef CONCEPT_SUBGRAPH_H
#define CONCEPT_SUBGRAPH_H


#include "Graph.h"


class Subgraph : public std::vector<Vertex> {
public:

    Subgraph(std::initializer_list<Vertex> list) : std::vector<Vertex>(list) {
#ifndef NDEBUG
        unsigned n = 0;
        for (auto x : list)
            for (auto y : list)
                n += (x == y);
        if (n != list.size()) std::cout << *this << "\n";
        assert(n == list.size() && "vertices are not unique");
#endif
    }

    class Vertices {
        const std::vector<Vertex> &m_vertices;
    public:
        explicit Vertices(const std::vector<Vertex> &vertices) : m_vertices(vertices) {}

        [[nodiscard]] auto begin() const { return m_vertices.begin(); }

        [[nodiscard]] auto end() const { return m_vertices.end(); }
    };

    [[nodiscard]] Vertices vertices() const {
        return Vertices(*this);
    }

    class VertexPairs {
        class Iterator {
            const std::vector<Vertex> &m_vertices;
            size_t i, j;
        public:
            Iterator(const std::vector<Vertex> &vertices, size_t i_start, size_t j_start) : m_vertices(vertices),
                                                                                            i(i_start), j(j_start) {}

            VertexPair operator*() const { return {m_vertices[i], m_vertices[j]}; }

            Iterator &operator++() {
                j++;
                if (j == m_vertices.size()) {
                    i++;
                    j = i + 1;
                }
                return *this;
            }

            bool operator==(const Iterator &other) const { return std::tie(i, j) == std::tie(other.i, other.j); }

            bool operator!=(const Iterator &other) const { return !(*this == other); }
        };

        const std::vector<Vertex> &m_vertices;
    public:
        explicit VertexPairs(const std::vector<Vertex> &vertices) : m_vertices(vertices) {}

        [[nodiscard]] Iterator begin() const { return Iterator(m_vertices, 0, 1); }

        [[nodiscard]] Iterator end() const {
            if (!m_vertices.empty()) {
                return Iterator(m_vertices, m_vertices.size() - 1, m_vertices.size());
            } else {
                return Iterator(m_vertices, 0, 1);
            }
        }
    };

    [[nodiscard]] VertexPairs vertexPairs() const {
        return VertexPairs(*this);
    }

    class UnmarkedVertexPairs {
        class Iterator {
            const std::vector<Vertex> &m_vertices;
            const VertexPairMap<bool> &m_marked;
            size_t i, j;
        public:
            Iterator(const std::vector<Vertex> &vertices, const VertexPairMap<bool> &marked, size_t i_start, size_t j_start)
                    : m_vertices(vertices), m_marked(marked), i(i_start), j(j_start) {
                const auto n = m_vertices.size();
                if (i < n - 1 && marked[{m_vertices[i], m_vertices[j]}]) ++(*this);
            }

            VertexPair operator*() const { return {m_vertices[i], m_vertices[j]}; }

            Iterator &operator++() {
                const auto n = m_vertices.size();
                do {
                    j++;
                    if (j == n) {
                        i++;
                        j = i + 1;
                    }
                } while (((i < n - 1) && m_marked[{m_vertices[i], m_vertices[j]}]));
                if (i == n - 1) { j = n; }
                return *this;
            }

            bool operator==(const Iterator &other) const { return std::tie(i, j) == std::tie(other.i, other.j); }

            bool operator!=(const Iterator &other) const { return !(*this == other); }
        };

        const std::vector<Vertex> &m_vertices;
        const VertexPairMap<bool> &m_marked;
    public:
        UnmarkedVertexPairs(const std::vector<Vertex> &vertices, const VertexPairMap<bool> &marked) :
        m_vertices(vertices), m_marked(marked) {}

        [[nodiscard]] auto begin() const { return Iterator(m_vertices, m_marked, 0, 1); }

        [[nodiscard]] auto end() const { return Iterator(m_vertices, m_marked, m_vertices.size() - 1, m_vertices.size()); }
    };

    [[nodiscard]] UnmarkedVertexPairs unmarkedVertexPairs(const VertexPairMap<bool> &marked) const {
        return UnmarkedVertexPairs(*this, marked);
    }


    template<typename VertexPairCallback>
    bool for_all_vertex_pairs(VertexPairCallback callback) const {
        for (size_t i = 0; i < size(); ++i) {
            for (size_t j = i + 1; j < size(); ++j) {
                if (callback(VertexPair((*this)[i], (*this)[j]))) return true;
            }
        }
        return false;
    }

    template<typename VertexPairCallback>
    bool for_all_unmarked_vertex_pairs(const VertexPairMap<bool> &marked, VertexPairCallback callback) const {
        return for_all_vertex_pairs([&](VertexPair uv) {
            return !marked[uv] && callback(uv);
            // if (!marked[uv]) return callback(uv);
            //return false;
        });
    }

    friend std::ostream &operator<<(std::ostream &os, const Subgraph &subgraph) {
        os << "{";
        for (Vertex u : subgraph) os << " " << u;
        return os << " }";
    }
};

constexpr Cost invalid_cost = std::numeric_limits<Cost>::max();

Cost cost(const Subgraph &subgraph, const VertexPairMap<bool> &marked, const VertexPairMap<Cost> &costs) {
    Cost min_cost = invalid_cost;
    subgraph.for_all_unmarked_vertex_pairs(marked, [&](VertexPair uv) {
        min_cost = std::min(min_cost, costs[uv]);
        return false;
    });
    return min_cost;
}

void verify(const Subgraph &subgraph, const Graph &graph) {
    assert(subgraph.size() == 4);

    const Subgraph &S = subgraph;

    assert(S[0] != S[1]);
    assert(S[0] != S[2]);
    assert(S[0] != S[3]);
    assert(S[1] != S[2]);
    assert(S[1] != S[3]);
    assert(S[2] != S[3]);

    bool is_one = false;
    std::vector<size_t> map = {0, 1, 2, 3};
    do {
        size_t i = map[0], j = map[1], k = map[2], l = map[3];
        bool is_path =
                 graph.has_edge({S[i], S[j]}) &&
                !graph.has_edge({S[i], S[k]}) &&
                !graph.has_edge({S[i], S[l]}) &&
                 graph.has_edge({S[j], S[k]}) &&
                !graph.has_edge({S[j], S[l]}) &&
                 graph.has_edge({S[k], S[l]});
        bool is_cycle =
                 graph.has_edge({S[i], S[j]}) &&
                !graph.has_edge({S[i], S[k]}) &&
                 graph.has_edge({S[i], S[l]}) &&
                 graph.has_edge({S[j], S[k]}) &&
                !graph.has_edge({S[j], S[l]}) &&
                 graph.has_edge({S[k], S[l]});
        is_one |= is_path || is_cycle;
    } while (std::next_permutation(map.begin(), map.end()));

    if (!is_one) {
        std::cout << subgraph << " E = {";
        for (VertexPair uv : subgraph.vertexPairs()) {
            if (graph.has_edge(uv)) std::cout << " " << uv;
        }
        std::cout << " }\n";
    }
    assert(is_one);
}


#endif //CONCEPT_SUBGRAPH_H
