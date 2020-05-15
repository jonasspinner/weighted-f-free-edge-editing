#ifndef WEIGHTED_F_FREE_EDGE_EDITING__SUBGRAPH_H
#define WEIGHTED_F_FREE_EDGE_EDITING__SUBGRAPH_H

#include "robin_hood.h"

#include "Graph.h"
#include "VertexPairMap.h"


class Subgraph {
    std::vector<Vertex> m_vertices;
public:

    /**
     * An induced subgraph G[S]. A reference to the corresponding graph is not stored.
     *
     * The subgraph is identified by the vertices S \subseteq V that are part of it. The order of the vertices is preserved.
     *
     * @param list
     */
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

    /**
     * An optimized constructor for concatenating the vertices and the subgraph.
     *
     * Only allocates memory once.
     *
     * @param A
     * @param other
     * @param B
     */
    Subgraph(std::initializer_list<Vertex> A, const Subgraph &other, std::initializer_list<Vertex> B) {
        m_vertices.reserve(A.size() + other.size() + B.size());
        m_vertices.insert(m_vertices.end(), A);
        m_vertices.insert(m_vertices.end(), other.m_vertices.begin(), other.m_vertices.end());
        m_vertices.insert(m_vertices.end(), B);
    }

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

    /**
     * Returns the Vertex at the given index in the subgraph.
     *
     * @param index
     */
    [[nodiscard]] const Vertex &operator[](size_t index) const { return m_vertices[index]; }

    friend std::ostream &operator<<(std::ostream &os, const Subgraph &subgraph) {
        os << "{";
        for (Vertex u : subgraph.m_vertices) os << " " << u;
        os << " }";
        return os;
    }

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const Subgraph &subgraph) {
        using namespace YAML;
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

    [[nodiscard]] bool contains(Vertex u) const {
        return std::any_of(m_vertices.begin(), m_vertices.end(), [&](Vertex y) { return y == u; });
    }

    [[nodiscard]] bool contains(VertexPair uv) const {
        return contains(uv.u) && contains(uv.v);
    }

    friend struct std::hash<Subgraph>;
};

template <>
struct std::hash<Subgraph> {
    size_t operator()(const Subgraph &subgraph) const noexcept {
        // hash_bytes has `void const* ptr` as first parameter type.
        auto ptr = static_cast<void const*>(subgraph.m_vertices.data());
        auto len = subgraph.m_vertices.size() * sizeof(Vertex); // length of m_vertices in bytes.
        return robin_hood::hash_bytes(ptr, len);
    }
};

constexpr Cost invalid_cost = std::numeric_limits<Cost>::max();

Cost get_subgraph_cost(const Subgraph &subgraph, const VertexPairMap<bool> &marked, const VertexPairMap<Cost> &costs);



#endif //WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPH_H
