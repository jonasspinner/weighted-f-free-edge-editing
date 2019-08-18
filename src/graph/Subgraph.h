//
// Created by jonas on 03.07.19.
//

#ifndef CONCEPT_SUBGRAPH_H
#define CONCEPT_SUBGRAPH_H


#include "Graph.h"


class Subgraph {
    std::vector<Vertex> m_vertices;
    std::string m_tag;
public:

    Subgraph(std::initializer_list<Vertex> list) : m_vertices(list) {
#ifndef NDEBUG
        unsigned n = 0;
        for (auto x : list)
            for (auto y : list)
                n += (x == y);
        if (n != list.size()) std::cout << *this << "\n";
        assert(n == list.size() && "vertices are not unique");
#endif
    }

    Subgraph(std::vector<Vertex> &&vertices, std::string &&tag) : m_vertices(std::move(vertices)), m_tag(std::move(tag)) {}

    class Vertices {
        const std::vector<Vertex> &m_vertices;
    public:
        explicit Vertices(const std::vector<Vertex> &vertices) : m_vertices(vertices) {}

        [[nodiscard]] auto begin() const { return m_vertices.begin(); }

        [[nodiscard]] auto end() const { return m_vertices.end(); }
    };

    /**
     * Returns Iterator over all vertices.
     *
     * @return
     */
    [[nodiscard]] Vertices vertices() const {
        return Vertices(m_vertices);
    }

    class VertexPairs {
        class Iterator {
            //const std::vector<Vertex> &m_vertices;
            //size_t i, j;
            using VertexIt = std::vector<Vertex>::const_iterator;
            VertexIt uit, vit, end;
        public:
            using value_type = VertexPair;
            using difference_type = std::ptrdiff_t;
            using pointer = const VertexPair*;
            using reference = const VertexPair&;
            using iterator_category = std::forward_iterator_tag;

            Iterator(VertexIt begin, VertexIt end_) : uit(begin), vit(begin), end(end_) {
                if (vit != end) ++vit;
                if (vit == end) uit = end;
            }

            VertexPair operator*() const {
                assert(uit != end);
                assert(vit != end);
                return {*uit, *vit};
            }

            Iterator &operator++() {
                assert(uit != end);
                assert(vit != end);
                ++vit;
                if (vit == end) {
                    ++uit;
                    vit = uit;
                    if (vit != end) ++vit;
                }
                if (vit == end) uit = end;
                return *this;
            }

            bool operator==(const Iterator &other) const { return std::tie(uit, vit) == std::tie(other.uit, other.vit); }

            bool operator!=(const Iterator &other) const { return !(*this == other); }
        };

        const std::vector<Vertex> &m_vertices;
    public:
        explicit VertexPairs(const std::vector<Vertex> &vertices) : m_vertices(vertices) {}

        [[nodiscard]] Iterator begin() const { return Iterator(m_vertices.begin(), m_vertices.end()); }

        [[nodiscard]] Iterator end() const { return Iterator(m_vertices.end(), m_vertices.end()); }
    };

    /**
     * Returns Iterator over all pairs of vertices {u, v}.
     *
     * @return
     */
    [[nodiscard]] VertexPairs vertexPairs() const {
        return VertexPairs(m_vertices);
    }

    [[nodiscard]] size_t size() const { return m_vertices.size(); }
    const Vertex &operator[](size_t index) const { return m_vertices[index]; }


    template<typename VertexPairCallback>
    bool for_all_vertex_pairs(VertexPairCallback callback) const {
        for (size_t i = 0; i < m_vertices.size(); ++i) {
            for (size_t j = i + 1; j < m_vertices.size(); ++j) {
                if (callback(VertexPair(m_vertices[i], m_vertices[j]))) return true;
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
        for (Vertex u : subgraph.m_vertices) os << " " << u;
        os << " }";
        if (!subgraph.m_tag.empty()) os << "#" << subgraph.m_tag;
        return os;
    }

    void sort_vertices() {
        std::sort(m_vertices.begin(), m_vertices.end());
    }

    void push_back(Vertex u) {
        m_vertices.push_back(u);
    }

    void append(const Subgraph &other) {
        m_vertices.insert(m_vertices.end(), other.m_vertices.begin(), other.m_vertices.end());
    }

    bool operator==(const Subgraph &other) const {
        return m_vertices == other.m_vertices;
    }

    bool operator<(const Subgraph &other) const {
        return std::lexicographical_compare(m_vertices.begin(), m_vertices.end(), other.m_vertices.begin(), other.m_vertices.end());
    }

    friend void verify(const Subgraph &, const Graph&);
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
    assert(subgraph.m_vertices.size() == 4);

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
