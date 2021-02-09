
#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHP3_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHP3_H


#include "robin_hood.h"

#include "Subgraph.h"


namespace subgraph_iterators {
    class CenterP3Finder;
}

template<>
class SubgraphT<Options::FSG::P3> {
    using Vertices = std::array<Vertex, 3>;
    Vertices m_vertices;

    constexpr SubgraphT(Vertex a, Vertex b, Vertex c) noexcept : m_vertices({a, b, c}) {}
    constexpr SubgraphT(const Vertices &vertices) noexcept : m_vertices(vertices) {}
public:
    friend struct std::hash<SubgraphT<Options::FSG::P3>>;
    friend class CenterP3Finder;

    using Finder = subgraph_iterators::CenterP3Finder;

    static constexpr SubgraphT P3(const Vertices &vertices) noexcept {
        return SubgraphT{vertices};
    }

    [[nodiscard]] constexpr auto size() const noexcept {
        return m_vertices.size();
    }

    [[nodiscard]] constexpr const auto &vertices() const noexcept {
        return m_vertices;
    }

    class VertexPairIt {
        std::array<VertexPair, 3> m_pairs{};
    public:
        constexpr VertexPairIt(const Vertices &v) noexcept: m_pairs{VertexPair{v[0], v[1]}, {v[0], v[2]}, {v[1], v[2]}} {}
        [[nodiscard]] constexpr auto begin() const noexcept {
            return m_pairs.begin();
        }
        [[nodiscard]] constexpr auto end() const noexcept {
            return m_pairs.end();
        }
    };
    [[nodiscard]] constexpr auto vertex_pairs() const noexcept {
        return VertexPairIt{m_vertices};
    }

    [[nodiscard]] constexpr auto non_converting_edits() const noexcept {
        return vertex_pairs();
    }

    [[nodiscard]] Cost calculate_min_cost(const VertexPairMap<Cost> &costs, const VertexPairMap<bool> &marked) const noexcept {
        const auto &[a, b, c] = m_vertices;
        auto x = [&](VertexPair uv) -> Cost { return marked[uv] ? invalid_cost : costs[uv]; };
        return std::min({x({a, b}), x({a, c}), x({b, c})});
    }

    [[nodiscard]] bool operator==(const SubgraphT &other) const noexcept {
        return m_vertices == other.m_vertices;
    }

    [[nodiscard]] bool operator!=(const SubgraphT &other) const noexcept {
        return !(*this == other);
    }

    [[nodiscard]] bool operator<(const SubgraphT &other) const noexcept {
        return m_vertices < other.m_vertices;
    }

    [[nodiscard]] bool contains(Vertex u) const noexcept {
        return std::any_of(m_vertices.begin(), m_vertices.end(), [&](Vertex y) { return y == u; });
    }

    [[nodiscard]] bool contains(VertexPair uv) const noexcept {
        return contains(uv.u) && contains(uv.v);
    }

    friend std::ostream &operator<<(std::ostream &os, const SubgraphT &subgraph) {
        os << "P3{";
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

    [[nodiscard]] constexpr Vertex &operator[](std::size_t index) noexcept {
        assert(index < size());
        return m_vertices[index];
    }

    [[nodiscard]] constexpr Vertex operator[](std::size_t index) const noexcept {
        assert(index < size());
        return m_vertices[index];
    }

    static bool is_valid_P3(const Graph& graph, const Vertices &vertices) noexcept {
        auto [a, b, c] = vertices;
        auto e = [&](VertexPair uv) { return graph.has_edge(uv); };
        return e({a, b}) && !e({a, c}) &&  e({b, c});
    }

    static bool is_valid_P3(const Graph& graph, const Graph &forbidden_graph, const Vertices &vertices) noexcept {
        auto [a, b, c] = vertices;
        auto e = [&](VertexPair uv) { return graph.has_edge(uv); };
        auto f = [&](VertexPair uv) { return forbidden_graph.has_edge(uv); };
        return  e({a, b}) && !e({a, c}) &&  e({b, c})
            && !f({a, b}) && !f({a, c}) && !f({b, c});
    }
};

template <>
struct std::hash<SubgraphT<Options::FSG::P3>> {
    size_t operator()(const SubgraphT<Options::FSG::P3> &subgraph) const noexcept {
        // hash_bytes has `void const* ptr` as first parameter type.
        const auto ptr = static_cast<void const*>(subgraph.m_vertices.data());
        const auto len = subgraph.m_vertices.size() * sizeof(Vertex); // length of m_vertices in bytes.
        return robin_hood::hash_bytes(ptr, len);
    }
};

namespace subgraph_iterators {

class CenterP3Finder {
    /**
     * Only allocates if temporary adjacency array rows are too small for the current graph.
     *
     * The invariants for all subgraphs {a, b, c} are
     *  + P3: a < c
     */
    Graph::AdjRow A;
    Graph::AdjRow B;
    Graph::AdjRow C;

    static inline void init(Graph::AdjRow &row, Vertex neighbor, Vertex non_neighbor, const Graph& graph) {
        row = graph.adj(neighbor);
        row -= graph.adj(non_neighbor);
        row[non_neighbor] = false;
    }

    /**
     * For neighbor u and non neighbor v find x, such that ux is an edge, vx is a non-edge and both are not forbidden.
     * Assumes that uv is an edge and not forbidden.
     */
    static inline void init(Graph::AdjRow &row, Vertex neighbor, Vertex non_neighbor, const Graph& graph, const Graph &forbidden_graph) {
        row = graph.adj(neighbor);
        row -= graph.adj(non_neighbor);
        row -= forbidden_graph.adj(neighbor);
        row -= forbidden_graph.adj(non_neighbor);
        row[non_neighbor] = false;
    }
public:
    using Subgraph = SubgraphT<Options::FSG::P3>;

    template<class Callback>
    IterationExit find(const Graph &graph, Callback callback) {
        static_assert(std::is_invocable_r_v<IterationControl, Callback, const Subgraph &>,
                "Callback must have IterationControl(const Subgraph &) signature.");

        /** P_3: <a, b, c> **/
        for (auto a : graph.vertices()) {
            for (auto b : graph.neighbors(a)) {

                // TODO: Optimize C to set the range 0..a to zero so that a < c for all c.
                init(C, b, a, graph);

                for (auto c : Graph::iterate(C)) {
                    assert(a != b); assert(a != c); assert(b != c);
                    assert(Subgraph::is_valid_P3(graph, {a, b, c}));

                    if (a < c) {
                        if (callback(Subgraph::P3({a, b, c})) == IterationControl::Break)
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

        /** P_3: <a, b, c> **/
        for (auto a : graph.vertices()) {
            for (auto b : graph.neighbors(a)) {
                if (forbidden_graph.has_edge({a, b}))
                    continue;

                // TODO: Optimize C to set the range 0..a to zero so that a < c for all c.
                init(C, b, a, graph, forbidden_graph);

                for (auto c : Graph::iterate(C)) {
                    assert(a != b); assert(a != c); assert(b != c);
                    assert(Subgraph::is_valid_P3(graph, {a, b, c}));

                    if (a < c) {
                        if (callback(Subgraph::P3({a, b, c})) == IterationControl::Break)
                            return IterationExit::Break;
                    }
                }
            }
        }
        return IterationExit::Normal;
    }

    template<class Callback>
    IterationExit find_unique(const Graph &graph, Callback callback) {
        return find(graph, callback);
    }

    template<class Callback>
    IterationExit find_near(VertexPair uv, const Graph &graph, const Graph &forbidden_graph, Callback callback) {
        static_assert(std::is_invocable_r_v<IterationControl, Callback, const Subgraph &>,
                "Callback must have IterationControl(const Subgraph &) signature.");
        auto [u, v] = uv;

        auto ensure_direction = [](auto &subgraph) {
            if (subgraph[0] > subgraph[2]) {
                std::swap(subgraph[0], subgraph[2]);
            }
        };

        if (forbidden_graph.has_edge(uv))
            return IterationExit::Normal;

        if (graph.has_edge(uv)) {
            init(A, u, v, graph, forbidden_graph);
            for (auto a : Graph::iterate(A)) {
                std::array<Vertex, 3> vertices{a, u, v};
                ensure_direction(vertices);
                assert(Subgraph::is_valid_P3(graph, vertices));
                if (callback(Subgraph::P3(vertices)) == IterationControl::Break)
                    return IterationExit::Break;
            }

            init(B, v, u, graph, forbidden_graph);
            for (auto b : Graph::iterate(B)) {
                std::array<Vertex, 3> vertices{u, v, b};
                ensure_direction(vertices);
                assert(Subgraph::is_valid_P3(graph, vertices));
                if (callback(Subgraph::P3(vertices)) == IterationControl::Break)
                    return IterationExit::Break;
            }
        } else {
            A = graph.adj(u);
            A &= graph.adj(v);
            A -= forbidden_graph.adj(u);
            A -= forbidden_graph.adj(v);

            for (auto a : Graph::iterate(A)) {
                assert(Subgraph::is_valid_P3(graph, {u, a, v}));
                if (callback(Subgraph::P3({u, a, v})) == IterationControl::Break)
                    return IterationExit::Break;
            }
        }
        return IterationExit::Normal;
    }
};

}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHP3_H
