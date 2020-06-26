#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHC4P4_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHC4P4_H

#include "robin_hood.h"

#include "Subgraph.h"
#include "../finder/CenterC4P4.h"

class CenterC4P4Finder;

template<>
class SubgraphT<Options::FSG::C4P4> {
    enum class Type {
        C4, P4
    } m_type;
    using Vertices = std::array<Vertex, 4>;
    Vertices m_vertices;

    constexpr SubgraphT(Type type, Vertex a, Vertex b, Vertex c, Vertex d) noexcept : m_type(type), m_vertices({a, b, c, d}) {}
    constexpr SubgraphT(Type type, const Vertices &vertices) noexcept : m_type(type), m_vertices(vertices) {}
public:
    friend struct std::hash<SubgraphT<Options::FSG::C4P4>>;
    friend class CenterC4P4Finder;

    using Subgraph = SubgraphT<Options::FSG::C4P4>;
    using Finder = CenterC4P4Finder;

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
    [[nodiscard]] constexpr auto vertices() const noexcept {
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

            constexpr Iterator() : m_u(), m_v(), m_end() {}

            constexpr explicit Iterator(Vertices::const_iterator u, Vertices::const_iterator v, Vertices::const_iterator end) noexcept : m_u(u), m_v(v), m_end(end) {}

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
    [[nodiscard]] constexpr auto vertex_pairs() const noexcept {
        return VertexPairIt{m_vertices};
    }

    class NonConvertingEdits {
        const Vertices &m_vertices;
    public:
        class Iterator {
            Vertices::const_iterator m_u;
            Vertices::const_iterator m_v;
            Vertices::const_iterator m_begin;
            Vertices::const_iterator m_end;
        public:
            using value_type = VertexPair;
            using difference_type = std::ptrdiff_t;
            using pointer = const VertexPair *;
            using reference = const VertexPair &;
            using iterator_category = std::forward_iterator_tag;

            constexpr Iterator() : m_u(), m_v(), m_begin(), m_end() {};

            constexpr explicit Iterator(Vertices::const_iterator u, Vertices::const_iterator v,
                    Vertices::const_iterator begin, Vertices::const_iterator end) : m_u(u), m_v(v), m_begin(begin), m_end(end) {}

            constexpr VertexPair operator*() const { return {*m_u, *m_v}; }
            constexpr Iterator &operator++() {
                ++m_v;
                if (m_v == m_end) {
                    ++m_u;
                    m_v = m_u + 1;
                } else if (m_u == m_begin && m_v + 1 == m_end) {
                    ++m_u;
                    m_v = m_u + 1;
                }
                return *this;
            }

            [[nodiscard]] constexpr bool operator==(const Iterator &other) const { return m_u == other.m_u && m_v == other.m_v; }

            [[nodiscard]] constexpr bool operator!=(const Iterator &other) const { return !(*this == other); }
        };
        constexpr NonConvertingEdits(const Vertices &vertices) noexcept : m_vertices(vertices) {}
        [[nodiscard]] constexpr auto begin() const {
            auto begin = m_vertices.begin();
            return Iterator(begin, begin + 1, begin,m_vertices.end());
        }
        [[nodiscard]] constexpr auto end() const {
            auto end = m_vertices.end();
            return Iterator(end - 1, end, m_vertices.begin(), end);
        }
    };
    [[nodiscard]] constexpr auto non_converting_edits() const noexcept {
        return NonConvertingEdits{m_vertices};
    }

    [[nodiscard]] Cost calculate_min_cost(const VertexPairMap<Cost> &costs, const VertexPairMap<bool> &marked) const {
        Cost min_cost = invalid_cost;
        for (auto uv : non_converting_edits()) {
            if (!marked[uv]) {
                min_cost = std::min(min_cost, costs[uv]);
            }
        }
        return min_cost;
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

    [[nodiscard]] constexpr Vertex &operator[](std::size_t index) noexcept {
        assert(index < size());
        return m_vertices[index];
    }

    [[nodiscard]] constexpr Vertex operator[](std::size_t index) const noexcept {
        assert(index < size());
        return m_vertices[index];
    }

    static bool is_valid_C4(const Graph& graph, const Graph &forbidden_graph, const Vertices &vertices) {
        auto [a, b, c, d] = vertices;
        auto e = [&](VertexPair uv) { return graph.hasEdge(uv); };
        auto f = [&](VertexPair uv) { return forbidden_graph.hasEdge(uv); };
        bool valid = true;
        valid &=  e({a, b}) && !e({a, c}) &&  e({a, d}) &&  e({b, c}) && !e({b, d}) &&  e({c, d});
        valid &= !f({a, b}) && !f({a, c}) &&               !f({b, c}) && !f({b, d}) && !f({c, d});
        return valid;
    }

    static bool is_valid_P4(const Graph& graph, const Graph &forbidden_graph, const Vertices &vertices) {
        auto [a, b, c, d] = vertices;
        auto e = [&](VertexPair uv) { return graph.hasEdge(uv); };
        auto f = [&](VertexPair uv) { return forbidden_graph.hasEdge(uv); };
        bool valid = true;
        valid &=  e({a, b}) && !e({a, c}) && !e({a, d}) &&  e({b, c}) && !e({b, d}) &&  e({c, d});
        valid &= !f({a, b}) && !f({a, c}) &&               !f({b, c}) && !f({b, d}) && !f({c, d});
        return valid;
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

    static inline void init(Graph::AdjRow &row, Vertex neighbor, Vertex non_neighbor, const Graph& graph) {
        row = graph.m_adj[neighbor];
        row -= graph.m_adj[non_neighbor];
        row[non_neighbor] = false;
    }

    /**
     * For neighbor u and non neighbor v find x, such that ux is an edge, vx is a non-edge and both are not forbidden.
     * Assumes that uv is an edge and not forbidden.
     */
    static inline void init(Graph::AdjRow &row, Vertex neighbor, Vertex non_neighbor, const Graph& graph, const Graph &forbidden_graph) {
        row = graph.m_adj[neighbor];
        row -= graph.m_adj[non_neighbor];
        row -= forbidden_graph.m_adj[neighbor];
        row -= forbidden_graph.m_adj[non_neighbor];
        row[non_neighbor] = false;
    }
public:
    using Subgraph = SubgraphT<Options::FSG::C4P4>;

    template<class Callback>
    bool find(const Graph &graph, Callback callback) {
        static_assert(std::is_invocable_r_v<bool, Callback, const Subgraph>, "Callback must have bool(Subgraph) signature.");

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
        static_assert(std::is_invocable_r_v<bool, Callback, const Subgraph>, "Callback must have bool(Subgraph) signature.");

        for (auto [u, v] : graph.edges()) {
            if (forbidden_graph.hasEdge({u, v}))
                continue;
            init(A, u, v, graph, forbidden_graph);
            init(B, v, u, graph, forbidden_graph);
            for (auto a : Graph::iterate(A)) {
                for (auto b : Graph::iterate(B)) {
                    // ab is allowed to be forbidden
                    if (graph.hasEdge({a, b})) {
                        assert(Subgraph::is_valid_C4(graph, forbidden_graph, {a, u, v, b}));
                        if (callback(Subgraph::C4({a, u, v, b})))
                            return true;
                    } else {
                        assert(Subgraph::is_valid_P4(graph, forbidden_graph, {a, u, v, b}));
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
        static_assert(std::is_invocable_r_v<bool, Callback, const Subgraph>, "Callback must have bool(Subgraph) signature.");
        auto [u, v] = uv;

        auto ensure_direction = [](auto &vertices) {
            if (vertices[1] > vertices[2]) {
                std::swap(vertices[0], vertices[3]);
                std::swap(vertices[1], vertices[2]);
            }
        };

        if (graph.hasEdge(uv) && !forbidden_graph.hasEdge(uv)) {

            init(A, u, v, graph, forbidden_graph);
            init(B, v, u, graph, forbidden_graph);
            for (auto a : Graph::iterate(A)) {
                for (auto b : Graph::iterate(B)) {
                    if (graph.hasEdge({a, b})) {
                        assert(Subgraph::is_valid_C4(graph, forbidden_graph, {a, u, v, b}));
                        if (callback(Subgraph::C4({a, u, v, b})))
                            return true;

                        if (!forbidden_graph.hasEdge({a, b})) {
                            Subgraph::Vertices vertices{u, a, b, v};
                            ensure_direction(vertices);
                            assert(Subgraph::is_valid_C4(graph, forbidden_graph, vertices));
                            if (callback(Subgraph::C4(vertices)))
                                return true;
                        }
                    } else {
                        assert(Subgraph::is_valid_P4(graph, forbidden_graph, {a, u, v, b}));
                        if (callback(Subgraph::P4({a, u, v, b})))
                            return true;
                    }
                }

                init(C, a, u, graph, forbidden_graph);
                for (auto c : Graph::iterate(C)) {
                    Subgraph::Vertices vertices{c, a, u, v};
                    ensure_direction(vertices);
                    if (graph.hasEdge({c, v})) {
                        assert(Subgraph::is_valid_C4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::C4(vertices)))
                            return true;
                    } else {
                        assert(Subgraph::is_valid_P4(graph, forbidden_graph, vertices));
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
                        assert(Subgraph::is_valid_C4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::C4(vertices)))
                            return true;
                    } else {
                        assert(Subgraph::is_valid_P4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::P4(vertices)))
                            return true;
                    }
                }
            }
        } else if (!graph.hasEdge(uv) && !forbidden_graph.hasEdge(uv)) {
            // a-u-b-v-c

            B = graph.m_adj[u];
            B &= graph.m_adj[v];
            B -= forbidden_graph.m_adj[u];
            B -= forbidden_graph.m_adj[v];
            B[u] = false;
            B[v] = false;


            A = graph.m_adj[u];
            A -= forbidden_graph.m_adj[u];

            C = graph.m_adj[v];
            C -= forbidden_graph.m_adj[v];

            for (auto b : Graph::iterate(B)) {
                for (auto a : Graph::iterate(A)) {
                    if (a == b || graph.hasEdge({a, b}) || forbidden_graph.hasEdge({a, b}))
                        continue;
                    Subgraph::Vertices vertices{a, u, b, v};
                    ensure_direction(vertices);
                    if (graph.hasEdge({a, v})) {
                        assert(Subgraph::is_valid_C4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::C4(vertices)))
                            return true;
                    } else {
                        assert(Subgraph::is_valid_P4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::P4(vertices)))
                            return true;
                    }
                }
                for (auto c : Graph::iterate(C)) {
                    if (c == b || graph.hasEdge({b, c}) || forbidden_graph.hasEdge({b, c}))
                        continue;
                    Subgraph::Vertices vertices{u, b, v, c};
                    ensure_direction(vertices);
                    if (graph.hasEdge({u, c})) {
                        assert(Subgraph::is_valid_C4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::C4(vertices)))
                            return true;
                    } else {
                        assert(Subgraph::is_valid_P4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::P4(vertices)))
                            return true;
                    }
                }
            }

            // u-a-c-v
            A -= B;
            C -= B;
            for (auto a : Graph::iterate(A)) {
                for (auto c : Graph::iterate(C)) {
                    if (forbidden_graph.hasEdge({a, v}) || forbidden_graph.hasEdge({c, u}) || forbidden_graph.hasEdge({a, c}))
                        continue;
                    if (graph.hasEdge({a, c})) {
                        Subgraph::Vertices vertices{u, a, c, v};
                        ensure_direction(vertices);
                        if (callback(Subgraph::P4(vertices)))
                            return true;
                    }
                }
            }
        } else if (forbidden_graph.hasEdge(uv)) {
            init(A, u, v, graph, forbidden_graph);
            init(B, v, u, graph, forbidden_graph);

            for (auto a : Graph::iterate(A)) {
                for (auto b : Graph::iterate(B)) {
                    if (graph.hasEdge({a, b}) && !forbidden_graph.hasEdge({a, b})) {
                        Subgraph::Vertices vertices{u, a, b, v};
                        ensure_direction(vertices);
                        if (graph.hasEdge(uv)) {
                            assert(Subgraph::is_valid_C4(graph, forbidden_graph, vertices));
                            if (callback(Subgraph::C4(vertices)))
                                return true;
                        } else {
                            assert(Subgraph::is_valid_P4(graph, forbidden_graph, vertices));
                            if (callback(Subgraph::P4(vertices)))
                                return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    template<class Callback>
    bool find_unique(const Graph &graph, Callback callback) {
        static_assert(std::is_invocable_r_v<bool, Callback, const Subgraph>, "Callback must have bool(Subgraph) signature.");

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
