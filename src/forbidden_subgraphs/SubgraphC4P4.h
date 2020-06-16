#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHC4P4_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHC4P4_H

#include "robin_hood.h"

#include "Subgraph.h"

class CenterC4P4Finder;

template<>
class SubgraphT<Options::FSG::C4P4> {
    enum class Type {
        C4, P4
    } m_type;
    using Vertices = std::array<Vertex, 4>;
    Vertices m_vertices;
public:
    friend struct std::hash<SubgraphT<Options::FSG::C4P4>>;
    friend class CenterC4P4Finder;

    using Subgraph = SubgraphT<Options::FSG::C4P4>;
    using Finder = CenterC4P4Finder;

    constexpr SubgraphT(Type type, Vertex a, Vertex b, Vertex c, Vertex d) noexcept : m_type(type), m_vertices({a, b, c, d}) {}
    constexpr SubgraphT(Type type, const Vertices &vertices) noexcept : m_type(type), m_vertices(vertices) {}

    static constexpr SubgraphT C4(const Vertices &vertices) noexcept {
        return SubgraphT{Type::C4, vertices};
    }

    static constexpr SubgraphT P4(const Vertices &vertices) noexcept {
        return SubgraphT{Type::P4, vertices};
    }

    [[nodiscard]] constexpr auto size() const noexcept {
        return m_vertices.size();
    }

    class VertexIt {
        const Vertices &m_vertices;
    public:
        constexpr VertexIt(const Vertices &vertices) noexcept : m_vertices(vertices) {}
        [[nodiscard]] constexpr auto begin() const {
            return m_vertices.begin();
        }
        [[nodiscard]] constexpr auto end() const {
            return m_vertices.end();
        }
    };
    [[nodiscard]] constexpr VertexIt vertices() const noexcept {
        return VertexIt{m_vertices};
    }

    class VertexPairIt {
        const Vertices &m_vertices;
    public:
        class Iterator {
            Vertices::const_iterator m_u;
            Vertices::const_iterator m_v;
            Vertices::const_iterator m_end;
        public:
            using value_type = VertexPair;
            using difference_type = std::ptrdiff_t;
            using pointer = const VertexPair *;
            using reference = const VertexPair &;
            using iterator_category = std::forward_iterator_tag;

            constexpr explicit Iterator(Vertices::const_iterator u, Vertices::const_iterator v, Vertices::const_iterator end) : m_u(u), m_v(v), m_end(end) {}

            constexpr VertexPair operator*() const { return {*m_u, *m_v}; }
            constexpr Iterator &operator++() {
                ++m_v;
                if (m_v == m_end) {
                    ++m_u;
                    m_v = m_u + 1;
                }
                return *this;
            }

            [[nodiscard]] constexpr bool operator==(const Iterator &other) const { return m_u == other.m_u && m_v == other.m_v; }

            [[nodiscard]] constexpr bool operator!=(const Iterator &other) const { return !(*this == other); }
        };
        constexpr VertexPairIt(const Vertices &vertices) noexcept : m_vertices(vertices) {}
        [[nodiscard]] constexpr auto begin() const {
            auto begin = m_vertices.begin();
            return Iterator(begin, begin + 1, m_vertices.end());
        }
        [[nodiscard]] constexpr auto end() const {
            auto end = m_vertices.end();
            return Iterator(end - 1, end, end);
        }
    };
    [[nodiscard]] constexpr VertexPairIt vertex_pairs() const noexcept {
        return VertexPairIt{m_vertices};
    }

    class NonConvertingEdits {
        const Vertices &m_vertices;
    public:
        constexpr NonConvertingEdits(const Vertices &vertices) noexcept : m_vertices(vertices) {}
        [[nodiscard]] constexpr auto begin() const {
            auto begin = m_vertices.begin();
            return VertexPairIt::Iterator(begin, begin + 1, m_vertices.end());
        }
        [[nodiscard]] constexpr auto end() const {
            auto end = m_vertices.end();
            return VertexPairIt::Iterator(end - 2, end - 1, end);
        }
    };
    [[nodiscard]] constexpr NonConvertingEdits non_converting_edits() const noexcept {
        return NonConvertingEdits{m_vertices};
    }

    [[nodiscard]] bool operator==(const Subgraph &other) const noexcept {
        return m_vertices == other.m_vertices;
    }

    [[nodiscard]] bool operator!=(const Subgraph &other) const noexcept {
        return !(*this == other);
    }

    [[nodiscard]] bool operator<(const Subgraph &other) const noexcept {
        return m_vertices < other.m_vertices;
    }

    [[nodiscard]] bool contains(Vertex u) const noexcept {
        return std::any_of(m_vertices.begin(), m_vertices.end(), [&](Vertex y) { return y == u; });
    }

    [[nodiscard]] bool contains(VertexPair uv) const noexcept {
        return contains(uv.u) && contains(uv.v);
    }

    friend std::ostream &operator<<(std::ostream &os, const Subgraph &subgraph) {
        switch(subgraph.m_type) {
            case Type::C4:
                os << "C4{";
                break;
            case Type::P4:
                os << "P4{";
                break;
        }
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

private:
    [[nodiscard]] constexpr Vertex &operator[](std::size_t index) noexcept {
        assert(index < size());
        return m_vertices[index];
    }

    [[nodiscard]] constexpr Vertex operator[](std::size_t index) const noexcept {
        assert(index < size());
        return m_vertices[index];
    }
};

template <>
struct std::hash<SubgraphT<Options::FSG::C4P4>> {
    size_t operator()(const SubgraphT<Options::FSG::C4P4> &subgraph) const noexcept {
        // hash_bytes has `void const* ptr` as first parameter type.
        auto ptr = static_cast<void const*>(subgraph.m_vertices.data());
        auto len = subgraph.m_vertices.size() * sizeof(Vertex); // length of m_vertices in bytes.
        return robin_hood::hash_bytes(ptr, len);
    }
};


class CenterC4P4Finder {
    /**
     * Only allocates if temporary adjacency array rows are too small for the current graph.
     *
     * The invariants for all subgraphs {a, b, c, d} are
     *  + for C4, P4: b < c and
     *  + for C4: b is minimum, a < c.
     */
    Graph::AdjRow A;
    Graph::AdjRow B;
    Graph::AdjRow C;

    inline void init(Graph::AdjRow &row, Vertex neighbor, Vertex non_neighbor, const Graph& graph) const {
        row.reset();
        row |= graph.m_adj[neighbor];
        row -= graph.m_adj[non_neighbor];
        row[non_neighbor] = false;
    }

    /**
     * For neighbor u and non neighbor v find x, such that ux is an edge, vx is a non-edge and both are not forbidden.
     * Assumes that uv is an edge and not forbidden.
     */
    inline void init(Graph::AdjRow &row, Vertex neighbor, Vertex non_neighbor, const Graph& graph, const Graph &forbidden_graph) const {
        row.reset();
        row |= graph.m_adj[neighbor];
        row -= graph.m_adj[non_neighbor];
        row -= forbidden_graph.m_adj[neighbor];
        row -= forbidden_graph.m_adj[non_neighbor];
        row[non_neighbor] = false;
    }
public:
    using Subgraph = SubgraphT<Options::FSG::C4P4>;

    template<class Callback>
    bool find(const Graph &graph, Callback callback) {
        static_assert(std::is_invocable_r_v<bool, Callback, Subgraph>, "Callback must have bool(Subgraph) signature.");

        A.resize(graph.size());
        B.resize(graph.size());

        for (auto [u, v] : graph.edges()) {
            init(A, u, v, graph);
            init(B, v, u, graph);
            for (auto a : Graph::iterate(A)) {
                for (auto b : Graph::iterate(B)) {
                    if (graph.hasEdge({a, b})) {
                        if (callback(Subgraph::C4({a, u, v, b})))
                            return true;
                    } else {
                        if (callback(Subgraph::P4({a, u, v, b})))
                            return true;
                    }
                }
            }
        }
        return false;
    }

    template<class Callback>
    bool find(const Graph &graph, const Graph &forbidden_graph, Callback callback) {
        static_assert(std::is_invocable_r_v<bool, Callback, Subgraph>, "Callback must have bool(Subgraph) signature.");

        A.resize(graph.size());
        B.resize(graph.size());

        for (auto [u, v] : graph.edges()) {
            if (forbidden_graph.hasEdge({u, v}))
                continue;
            init(A, u, v, graph, forbidden_graph);
            init(B, v, u, graph, forbidden_graph);
            for (auto a : Graph::iterate(A)) {
                for (auto b : Graph::iterate(B)) {
                    // ab is allowed to be forbidden
                    if (graph.hasEdge({a, b})) {
                        if (callback(Subgraph::C4({a, u, v, b})))
                            return true;
                    } else {
                        if (callback(Subgraph::P4({a, u, v, b})))
                            return true;
                    }
                }
            }
        }
        return false;
    }

    template<class Callback>
    bool find_near(VertexPair uv, const Graph &graph, const Graph &forbidden_graph, Callback callback) {
        static_assert(std::is_invocable_r_v<bool, Callback, Subgraph>, "Callback must have bool(Subgraph) signature.");
        auto [u, v] = uv;

        A.resize(graph.size());
        B.resize(graph.size());
        C.resize(graph.size());

        auto ensure_direction = [](auto &vertices) {
            if (vertices[1] > vertices[2]) {
                std::swap(vertices[0], vertices[3]);
                std::swap(vertices[1], vertices[2]);
            }
        };


        if (forbidden_graph.hasEdge({u, v}))
            return false;

        init(A, u, v, graph, forbidden_graph);
        init(B, v, u, graph, forbidden_graph);
        for (auto a : Graph::iterate(A)) {
            for (auto b : Graph::iterate(B)) {
                if (graph.hasEdge({a, b})) {
                    if (callback(Subgraph::C4({a, u, v, b})))
                        return true;
                } else {
                    if (callback(Subgraph::P4({a, u, v, b})))
                        return true;
                }
            }

            init(C, a, u, graph, forbidden_graph);
            for (auto c : Graph::iterate(C)) {
                Subgraph::Vertices vertices{c, a, u, v};
                ensure_direction(vertices);
                if (graph.hasEdge({c, v})) {
                    if (callback(Subgraph::C4(vertices)))
                        return true;
                } else {
                    if (callback(Subgraph::P4(vertices)))
                        return true;
                }
            }
        }

        for (auto b : Graph::iterate(B)) {
            init(C, b, v, graph, forbidden_graph);
            for (auto c : Graph::iterate(C)) {
                Subgraph::Vertices vertices{u, v, b, c};
                ensure_direction(vertices);
                if (graph.hasEdge({c, u})) {
                    if (callback(Subgraph::C4(vertices)))
                        return true;
                } else {
                    if (callback(Subgraph::P4(vertices)))
                        return true;
                }
            }
        }
        return false;
    }

    template<class Callback>
    bool find_unique(const Graph &graph, Callback callback) {
        static_assert(std::is_invocable_r_v<bool, Callback, Subgraph>, "Callback must have bool(Subgraph) signature.");

        auto is_correct_cycle = [](const auto &subgraph) {
            const auto [a, b, c, d] = subgraph.m_vertices;
            return b < std::min({a, c, d}) && a < c;
        };

        return find(graph, [&](Subgraph subgraph) {
            if (subgraph.m_type == Subgraph::Type::C4) {
                if (!is_correct_cycle(subgraph))
                    return false;
            }
            if (callback(subgraph))
                return true;
            return false;
        });
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHC4P4_H
