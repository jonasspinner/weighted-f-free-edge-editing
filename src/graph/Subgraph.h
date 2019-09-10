//
// Created by jonas on 03.07.19.
//

#ifndef CONCEPT_SUBGRAPH_H
#define CONCEPT_SUBGRAPH_H


#include "Graph.h"
#include "VertexPairMap.h"


class Subgraph {
    std::vector<Vertex> m_vertices;
#ifndef NDEBUG
    std::string m_tag;
#endif
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

    Subgraph(std::initializer_list<Vertex> A, const Subgraph &other, std::initializer_list<Vertex> B) {
        m_vertices.reserve(A.size() + other.size() + B.size());
        m_vertices.insert(m_vertices.end(), A);
        m_vertices.insert(m_vertices.end(), other.m_vertices.begin(), other.m_vertices.end());
        m_vertices.insert(m_vertices.end(), B);
    }

#ifndef NDEBUG
    Subgraph(std::vector<Vertex> &&vertices, std::string &&tag) : m_vertices(std::move(vertices)), m_tag(std::move(tag)) {}
#else
    [[deprecated]] Subgraph(std::vector<Vertex> &&vertices, std::string &&tag) : m_vertices(std::move(vertices)) {}
#endif

    class Vertices {
        const std::vector<Vertex> &m_vertices;
    public:
        explicit Vertices(const std::vector<Vertex> &vertices) : m_vertices(vertices) {}

        [[nodiscard]] auto begin() const { return m_vertices.begin(); }

        [[nodiscard]] auto end() const { return m_vertices.end(); }
    };

    /**
     * Returns a range over the vertices in the subgraph.
     *
     * @return
     */
    [[nodiscard]] Vertices vertices() const {
        return Vertices(m_vertices);
    }

    class VertexPairs {
        class Iterator {
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
     * Returns a range over all pairs of vertices {u, v} in the subgraph.
     *
     * @return
     */
    [[nodiscard]] VertexPairs vertexPairs() const {
        return VertexPairs(m_vertices);
    }

    [[nodiscard]] size_t size() const { return m_vertices.size(); }
    const Vertex &operator[](size_t index) const { return m_vertices[index]; }

    friend std::ostream &operator<<(std::ostream &os, const Subgraph &subgraph) {
        os << "{";
        for (Vertex u : subgraph.m_vertices) os << " " << u;
        os << " }";
#ifndef NDEBUG
        if (!subgraph.m_tag.empty()) os << "#" << subgraph.m_tag;
#endif
        return os;
    }

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const Subgraph &subgraph) {
        out << YAML::Flow << YAML::BeginSeq;
        for (Vertex u : subgraph.vertices()) out << u;
        return out << YAML::EndSeq;
    }

    void sortVertices() {
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

Cost get_subgraph_cost(const Subgraph &subgraph, const VertexPairMap<bool> &marked, const VertexPairMap<Cost> &costs);

void verify(const Subgraph &subgraph, const Graph &graph);


#endif //CONCEPT_SUBGRAPH_H
