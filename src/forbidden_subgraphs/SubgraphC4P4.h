#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHC4P4_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHC4P4_H

#include "robin_hood.h"

#include "Subgraph.h"
#include "IndexedVertexPairRange.h"


namespace subgraph_iterators {
    // Forward declaration.
    class CenterC4P4Finder;
}

template<>
class SubgraphT<Options::FSG::C4P4> {
    friend class subgraph_iterators::CenterC4P4Finder;

    friend struct std::hash<SubgraphT<Options::FSG::C4P4>>;

    enum class Type {
        C4, P4
    } m_type;
    using Vertices = std::array<Vertex, 4>;
    Vertices m_vertices;

    SubgraphT() noexcept = default;

    constexpr SubgraphT(Type type, const Vertices &vertices): m_type(type), m_vertices(vertices) {
        assert(m_vertices[1] < m_vertices[2]);
    }

    constexpr SubgraphT(Type type, Vertices &&vertices): m_type(type), m_vertices(std::move(vertices)) {
        assert(m_vertices[1] < m_vertices[2]);
    }

public:
    using Finder = subgraph_iterators::CenterC4P4Finder;

    static constexpr auto C4(const Vertices &vertices) {
        return SubgraphT{Type::C4, vertices};
    }

    static constexpr auto C4(Vertices &&vertices) {
        return SubgraphT{Type::C4, std::move(vertices)};
    }

    static constexpr auto P4(const Vertices &vertices) {
        return SubgraphT{Type::P4, vertices};
    }

    static constexpr auto P4(Vertices &&vertices) {
        return SubgraphT{Type::P4, std::move(vertices)};
    }

    void swap(SubgraphT &other) noexcept {
        using std::swap;
        swap(m_type, other.m_type);
        swap(m_vertices, other.m_vertices);
    }

    [[nodiscard]] constexpr auto size() const noexcept {
        return m_vertices.size();
    }

    [[nodiscard]] constexpr const auto &vertices() const noexcept {
        return m_vertices;
    }

private:
    static constexpr std::array<IndexedVertexPairRange::IndexPair, 6> m_pair_indices{std::pair{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
public:
    [[nodiscard]] constexpr auto vertex_pairs() const noexcept {
        return IndexedVertexPairRange{m_vertices, m_pair_indices};
    }

private:
    static constexpr std::array<IndexedVertexPairRange::IndexPair, 5> m_non_converting_edits_indices{std::pair{0, 1}, {0, 2}, {1, 2}, {1, 3}, {2, 3}};
public:
    [[nodiscard]] constexpr auto non_converting_edits() const noexcept {
        return IndexedVertexPairRange{m_vertices, m_non_converting_edits_indices};
    }

    [[nodiscard]] Cost calculate_min_cost(const VertexPairMap<Cost> &costs, const VertexPairMap<bool> &marked) const {
        const auto &[a, b, c, d] = m_vertices;
        auto x = [&](VertexPair uv) -> Cost { return marked[uv] ? invalid_cost : costs[uv]; };
        return std::min({x({a, b}), x({a, c}), x({b, c}), x({b, d}), x({c, d})});
    }

    [[nodiscard]] constexpr bool operator==(const SubgraphT &other) const noexcept {
        return m_type == other.m_type && m_vertices == other.m_vertices;
    }

    [[nodiscard]] constexpr bool operator!=(const SubgraphT &other) const noexcept {
        return !(*this == other);
    }

    [[nodiscard]] bool operator<(const SubgraphT &other) const noexcept {
        return m_vertices < other.m_vertices;
    }

    [[nodiscard]] constexpr bool contains(Vertex u) const noexcept {
        const auto &[a, b, c, d] = m_vertices;
        return (a == u) || (b == u) || (c == u) || (d == u);
        // return std::any_of(m_vertices.begin(), m_vertices.end(), [&](Vertex y) { return y == u; }); // not constexpr
    }

    [[nodiscard]] constexpr bool contains(VertexPair uv) const noexcept {
        return contains(uv.u) && contains(uv.v);
    }

    friend std::ostream &operator<<(std::ostream &os, const SubgraphT &subgraph) {
        switch (subgraph.m_type) {
            case Type::C4:
                os << "C4{";
                break;
            case Type::P4:
                os << "P4{";
                break;
        }
        for (Vertex u : subgraph.vertices())
            os << " " << u;
        os << " }";
        return os;
    }

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const SubgraphT &subgraph) {
        using namespace YAML;
        out << YAML::Flow << YAML::BeginSeq;
        for (Vertex u : subgraph.vertices())
            out << u;
        return out << YAML::EndSeq;
    }

    [[nodiscard]] constexpr Vertex &operator[](std::size_t index) {
        assert(index < size());
        return m_vertices[index];
    }

    [[nodiscard]] constexpr Vertex operator[](std::size_t index) const {
        assert(index < size());
        return m_vertices[index];
    }

    static bool is_valid_C4(const Graph &graph, const Graph &forbidden_graph, const Vertices &vertices) {
        const auto[a, b, c, d] = vertices;
        const auto e = [&](VertexPair uv) { return graph.has_edge(uv); };
        const auto f = [&](VertexPair uv) { return forbidden_graph.has_edge(uv); };
        bool valid = true;
        valid &=  e({a, b}) && !e({a, c}) &&  e({a, d}) &&  e({b, c}) && !e({b, d}) &&  e({c, d});
        valid &= !f({a, b}) && !f({a, c}) &&               !f({b, c}) && !f({b, d}) && !f({c, d});
        return valid;
    }

    static bool is_valid_P4(const Graph &graph, const Graph &forbidden_graph, const Vertices &vertices) {
        const auto[a, b, c, d] = vertices;
        const auto e = [&](VertexPair uv) { return graph.has_edge(uv); };
        const auto f = [&](VertexPair uv) { return forbidden_graph.has_edge(uv); };
        bool valid = true;
        valid &=  e({a, b}) && !e({a, c}) && !e({a, d}) &&  e({b, c}) && !e({b, d}) &&  e({c, d});
        valid &= !f({a, b}) && !f({a, c}) &&               !f({b, c}) && !f({b, d}) && !f({c, d});
        return valid;
    }
};

template<>
struct std::hash<SubgraphT<Options::FSG::C4P4>> {
    inline size_t operator()(const SubgraphT<Options::FSG::C4P4> &subgraph) const noexcept {
        // hash_bytes has `void const* ptr` as first parameter type.
        const auto ptr = static_cast<void const *>(subgraph.m_vertices.data());
        const auto len = subgraph.m_vertices.size() * sizeof(Vertex); // length of m_vertices in bytes.
        return robin_hood::hash_bytes(ptr, len);
    }
};

template<>
inline void std::swap<SubgraphT<Options::FSG::C4P4>>(SubgraphT<Options::FSG::C4P4> &lhs, SubgraphT<Options::FSG::C4P4> &rhs) noexcept {
    lhs.swap(rhs);
}

namespace subgraph_iterators {

class CenterC4P4Finder {
    /**
     * find_near with forbidden graph: LocalSearch
     * find_near_unique without forbidden graph: SubgraphStats
     */
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

    static inline void init(Graph::AdjRow &row, Vertex neighbor, Vertex non_neighbor, const Graph &graph) {
        row = graph.adj(neighbor);
        row -= graph.adj(non_neighbor);
        row[non_neighbor] = false;
    }

    /**
     * For neighbor u and non neighbor v find x, such that ux is an edge, vx is a non-edge and both are not forbidden.
     * Assumes that uv is an edge and not forbidden.
     */
    static inline void
    init(Graph::AdjRow &row, Vertex neighbor, Vertex non_neighbor, const Graph &graph, const Graph &forbidden_graph) {
        row = graph.adj(neighbor);
        row -= graph.adj(non_neighbor);
        row -= forbidden_graph.adj(neighbor);
        row -= forbidden_graph.adj(non_neighbor);
        row[non_neighbor] = false;
    }

public:
    using Subgraph = SubgraphT<Options::FSG::C4P4>;

    template<class Callback>
    IterationExit find(const Graph &graph, Callback callback) {
        static_assert(std::is_invocable_r_v<IterationControl, Callback, const Subgraph &>,
                      "Callback must have IterationControl(const Subgraph &) signature.");

        for (const auto[u, v] : graph.edges()) {
            init(A, u, v, graph);
            init(B, v, u, graph);
            for (auto a : Graph::iterate(A)) {
                for (auto b : Graph::iterate(B)) {
                    if (graph.has_edge({a, b})) {
                        if (callback(Subgraph::C4({a, u, v, b})) == IterationControl::Break)
                            return IterationExit::Break;
                    } else {
                        if (callback(Subgraph::P4({a, u, v, b})) == IterationControl::Break)
                            return IterationExit::Break;
                    }
                }
            }
        }
        return IterationExit::Normal;
    }

    template<class Callback>
    IterationExit find(const Graph &graph, const Graph &forbidden_graph, Callback callback) {
        static_assert(std::is_invocable_r_v<IterationControl, Callback, const Subgraph &>,
                      "Callback must have IterationControl(const Subgraph &) signature.");

        for (auto uv : graph.edges()) {
            if (forbidden_graph.has_edge(uv))
                continue;
            const auto[u, v] = uv;
            init(A, u, v, graph, forbidden_graph);
            init(B, v, u, graph, forbidden_graph);
            for (auto a : Graph::iterate(A)) {
                for (auto b : Graph::iterate(B)) {
                    // ab is allowed to be forbidden
                    if (graph.has_edge({a, b})) {
                        assert(Subgraph::is_valid_C4(graph, forbidden_graph, {a, u, v, b}));
                        if (callback(Subgraph::C4({a, u, v, b})) == IterationControl::Break)
                            return IterationExit::Break;
                    } else {
                        assert(Subgraph::is_valid_P4(graph, forbidden_graph, {a, u, v, b}));
                        if (callback(Subgraph::P4({a, u, v, b})) == IterationControl::Break)
                            return IterationExit::Break;
                    }
                }
            }
        }
        return IterationExit::Normal;
    }

    template<class Callback>
    IterationExit find_near(VertexPair uv, const Graph &graph, const Graph &forbidden_graph, Callback callback) {
        static_assert(std::is_invocable_r_v<IterationControl, Callback, const Subgraph &>,
                      "Callback must have IterationControl(const Subgraph &) signature.");
        const auto[u, v] = uv;

        constexpr auto ensure_direction = [](auto &vertices) constexpr noexcept {
            if (vertices[1] > vertices[2]) {
                std::swap(vertices[0], vertices[3]);
                std::swap(vertices[1], vertices[2]);
            }
        };
#ifndef NDEBUG
        const auto e = [&](VertexPair pair) { return graph.has_edge(pair); };
        const auto f = [&](VertexPair pair) { return forbidden_graph.has_edge(pair); };
#endif

        if (graph.has_edge(uv) && !forbidden_graph.has_edge(uv)) {

            init(A, u, v, graph, forbidden_graph);
            init(B, v, u, graph, forbidden_graph);

            for (auto a : Graph::iterate(A)) {
                for (auto b : Graph::iterate(B)) {
                    if (graph.has_edge({a, b})) {
                        assert(Subgraph::is_valid_C4(graph, forbidden_graph, {a, u, v, b}));
                        if (callback(Subgraph::C4({a, u, v, b})) == IterationControl::Break)
                            return IterationExit::Break;

                        if (!forbidden_graph.has_edge({a, b})) {
                            Subgraph::Vertices vertices{u, a, b, v};
                            ensure_direction(vertices);
                            assert(Subgraph::is_valid_C4(graph, forbidden_graph, vertices));
                            if (callback(Subgraph::C4(vertices)) == IterationControl::Break)
                                return IterationExit::Break;
                        }
                    } else {
                        assert(Subgraph::is_valid_P4(graph, forbidden_graph, {a, u, v, b}));
                        if (callback(Subgraph::P4({a, u, v, b})) == IterationControl::Break)
                            return IterationExit::Break;
                    }
                }

                init(C, a, u, graph, forbidden_graph);
                for (auto c : Graph::iterate(C)) {
                    Subgraph::Vertices vertices{c, a, u, v};
                    ensure_direction(vertices);
                    if (graph.has_edge({c, v})) {
                        assert(Subgraph::is_valid_C4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::C4(vertices)) == IterationControl::Break)
                            return IterationExit::Break;
                    } else {
#ifndef NDEBUG
                        auto[v0, v1, v2, v3] = vertices;
                        auto e0 = e({v0, v1}), e1 = e({v0, v2}),   e2 = e({v0, v3}),   e3 = e({v1, v2}), e4 = e({v1, v3}), e5 = e({v2, v3});
                        auto f0 = f({v0, v1}), f1 = f({v0, v2}), /*f2 = f({v0, v3}),*/ f3 = f({v1, v2}), f4 = f({v1, v3}), f5 = f({v2, v3});
                        assert( e0); assert(!e1); assert(!e2); assert( e3); assert(!e4); assert( e5);
                        assert(!f0); assert(!f1);              assert(!f3); assert(!f4); assert(!f5);
#endif
                        assert(Subgraph::is_valid_P4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::P4(vertices)) == IterationControl::Break)
                            return IterationExit::Break;
                    }
                }
            }

            for (auto b : Graph::iterate(B)) {
                init(C, b, v, graph, forbidden_graph);
                for (auto c : Graph::iterate(C)) {
                    Subgraph::Vertices vertices{u, v, b, c};
                    ensure_direction(vertices);
                    if (graph.has_edge({c, u})) {
                        assert(Subgraph::is_valid_C4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::C4(vertices)) == IterationControl::Break)
                            return IterationExit::Break;
                    } else {
#ifndef NDEBUG
                        auto[v0, v1, v2, v3] = vertices;
                        auto e0 = e({v0, v1}), e1 = e({v0, v2}),   e2 = e({v0, v3}),   e3 = e({v1, v2}), e4 = e({v1, v3}), e5 = e({v2, v3});
                        auto f0 = f({v0, v1}), f1 = f({v0, v2}), /*f2 = f({v0, v3}),*/ f3 = f({v1, v2}), f4 = f({v1, v3}), f5 = f({v2, v3});
                        assert( e0); assert(!e1); assert(!e2); assert( e3); assert(!e4); assert( e5);
                        assert(!f0); assert(!f1);              assert(!f3); assert(!f4); assert(!f5);
#endif
                        assert(Subgraph::is_valid_P4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::P4(vertices)) == IterationControl::Break)
                            return IterationExit::Break;
                    }
                }
            }
        } else if (!graph.has_edge(uv) && !forbidden_graph.has_edge(uv)) {
            // a-u-b-v-c

            B = graph.adj(u);
            B &= graph.adj(v);
            B -= forbidden_graph.adj(u);
            B -= forbidden_graph.adj(v);
            B[u] = false;
            B[v] = false;


            A = graph.adj(u);
            A -= forbidden_graph.adj(u);

            C = graph.adj(v);
            C -= forbidden_graph.adj(v);

            for (auto b : Graph::iterate(B)) {
                for (auto a : Graph::iterate(A)) {
                    if (a == b || graph.has_edge({a, b}) || forbidden_graph.has_edge({a, b}))
                        continue;
                    Subgraph::Vertices vertices{a, u, b, v};
                    ensure_direction(vertices);
                    if (graph.has_edge({a, v})) {
                        assert(Subgraph::is_valid_C4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::C4(vertices)) == IterationControl::Break)
                            return IterationExit::Break;
                    } else {
                        assert(Subgraph::is_valid_P4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::P4(vertices)) == IterationControl::Break)
                            return IterationExit::Break;
                    }
                }
                for (auto c : Graph::iterate(C)) {
                    if (c == b || graph.has_edge({b, c}) || forbidden_graph.has_edge({b, c}))
                        continue;
                    Subgraph::Vertices vertices{u, b, v, c};
                    ensure_direction(vertices);
                    if (graph.has_edge({u, c})) {
                        assert(Subgraph::is_valid_C4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::C4(vertices)) == IterationControl::Break)
                            return IterationExit::Break;
                    } else {
                        assert(Subgraph::is_valid_P4(graph, forbidden_graph, vertices));
                        if (callback(Subgraph::P4(vertices)) == IterationControl::Break)
                            return IterationExit::Break;
                    }
                }
            }

            // u-a-c-v
            A -= B;
            A -= forbidden_graph.adj(v);
            C -= B;
            C -= forbidden_graph.adj(u);
            for (auto a : Graph::iterate(A)) {
                for (auto c : Graph::iterate(C)) {
                    if (forbidden_graph.has_edge({a, c}))
                        continue;
                    if (graph.has_edge({a, c})) {
                        Subgraph::Vertices vertices{u, a, c, v};
                        ensure_direction(vertices);
                        if (callback(Subgraph::P4(vertices)) == IterationControl::Break)
                            return IterationExit::Break;
                    }
                }
            }
        } else if (forbidden_graph.has_edge(uv)) {
            init(A, u, v, graph, forbidden_graph);
            init(B, v, u, graph, forbidden_graph);

            for (auto a : Graph::iterate(A)) {
                for (auto b : Graph::iterate(B)) {
                    if (graph.has_edge({a, b}) && !forbidden_graph.has_edge({a, b})) {
                        Subgraph::Vertices vertices{u, a, b, v};
                        ensure_direction(vertices);
                        if (graph.has_edge(uv)) {
                            assert(Subgraph::is_valid_C4(graph, forbidden_graph, vertices));
                            if (callback(Subgraph::C4(vertices)) == IterationControl::Break)
                                return IterationExit::Break;
                        } else {
                            assert(Subgraph::is_valid_P4(graph, forbidden_graph, vertices));
                            if (callback(Subgraph::P4(vertices)) == IterationControl::Break)
                                return IterationExit::Break;
                        }
                    }
                }
            }
        }
        return IterationExit::Normal;
    }

    template<class Callback>
    IterationExit find_unique(const Graph &graph, Callback callback) {
        static_assert(std::is_invocable_r_v<IterationControl, Callback, const Subgraph &>,
                      "Callback must have IterationControl(const Subgraph &) signature.");

        constexpr auto is_correct_cycle = [](const auto &subgraph) {
            const auto&[a, b, c, d] = subgraph.m_vertices;
            return b < std::min({a, c, d}) && a < c;
        };

        return find(graph, [&](const Subgraph &subgraph) {
            if (subgraph.m_type == Subgraph::Type::C4) {
                if (!is_correct_cycle(subgraph))
                    return IterationControl::Continue;
            }
            return callback(subgraph);
        });
    }

    template<class Callback>
    IterationExit find_near_unique(VertexPair uv, const Graph &graph, const Graph &forbidden_graph, Callback callback) {
        static_assert(std::is_invocable_r_v<IterationControl, Callback, const Subgraph>,
                      "Callback must have IterationControl(Subgraph) signature.");

        constexpr auto is_correct_cycle = [](const auto &subgraph) {
            const auto[a, b, c, d] = subgraph.m_vertices;
            return b < std::min({a, c, d}) && a < c;
        };

        return find_near(uv, graph, forbidden_graph, [&](const Subgraph &subgraph) {
            if (subgraph.m_type == Subgraph::Type::C4) {
                if (!is_correct_cycle(subgraph))
                    return IterationControl::Continue;
            }
            return callback(subgraph);
        });
    }
};

}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHC4P4_H
