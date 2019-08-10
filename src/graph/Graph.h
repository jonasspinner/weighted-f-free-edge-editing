//
// Created by jonas on 27.06.19.
//

#ifndef CONCEPT_GRAPH_H
#define CONCEPT_GRAPH_H

#include <iostream>
#include <vector>
#include <boost/dynamic_bitset.hpp>

#include "../definitions.h"

namespace Finder {
    class NaiveC4P4;

    class CenterC4P4;

    class CenterP3;

    class NaiveP3;
}
// template <int length> class Center;
namespace detail {
    template<int length>
    class Center;

    template<int length, bool with_cycles>
    class CenterFinderImpl;
}


class VertexPair {
public:
    Vertex u;
    Vertex v;

    VertexPair(Vertex x, Vertex y) : u(x < y ? x : y), v(x < y ? y : x) { assert(x != y); }

    friend std::ostream &operator<<(std::ostream &os, VertexPair uv) {
        return os << "{" << uv.u << ", " << uv.v << "}";
    }

    friend bool operator==(VertexPair uv, VertexPair xy) {
        return (uv.u == xy.u) && (uv.v == xy.v);
    }
};


class Graph {
public:
    using Block = uint_fast32_t;
    using AdjRow = boost::dynamic_bitset<Block>;
    using AdjMatrix = std::vector<AdjRow>;

private:
    unsigned int n;
    AdjMatrix adj;

public:

    explicit Graph(unsigned int size) : n(size), adj(n, boost::dynamic_bitset<Block>(n)) {}

    /**
     * Returns the number of vertices.
     *
     * @return
     */
    [[nodiscard]] Vertex size() const { return n; }


    void clear_edges() {
        for (auto &row: adj) { row.reset(); }
    }

    /**
     * Toggles the edge.
     *
     * @param edge
     */
    void toggle_edge(VertexPair edge) {
        const auto&[u, v] = edge;
        adj[u].flip(v);
        adj[v].flip(u);
    }

    /**
     * Returns the degree of vertex u.
     *
     * @param u
     * @return
     */
    [[nodiscard]] size_t degree(Vertex u) const { return adj[u].count(); };

    /**
     * Checks whether the edge is in the graph.
     *
     * @param edge
     * @return
     */
    [[nodiscard]] bool has_edge(VertexPair edge) const { return adj[edge.u][edge.v]; }

    /**
     * Sets the edge.
     *
     * @param edge
     */
    void set_edge(VertexPair edge) {
        const auto&[u, v] = edge;
        adj[u].set(v);
        adj[v].set(u);
    }

    void set_edges(const std::vector<VertexPair> &edges) {
        for (VertexPair uv : edges) {
            assert(!has_edge(uv));
            set_edge(uv);
        }
    }

    /**
     * Clears the edge.
     *
     * @param edge
     */
    void clear_edge(VertexPair edge) {
        const auto&[u, v] = edge;
        adj[u].reset(v);
        adj[v].reset(u);
    }


    class Vertices {
        class Iterator {
            Vertex v;
        public:
            explicit Iterator(Vertex start) : v(start) {}

            Vertex operator*() const { return v; }

            Iterator &operator++() {
                ++v;
                return *this;
            }

            bool operator==(const Iterator &other) const { return v == other.v; }

            bool operator!=(const Iterator &other) const { return !(*this == other); }
        };

        Vertex n;
    public:
        explicit Vertices(Vertex size) : n(size) {}

        [[nodiscard]] Iterator begin() const { return Iterator(0); }

        [[nodiscard]] Iterator end() const { return Iterator(n); }
    };

    [[nodiscard]] Vertices vertices() const {
        return Vertices(n);
    }


    class RowVertices {
        class Iterator {
            const AdjRow &m_row;
            Vertex u;
        public:
            explicit Iterator(const AdjRow &row) : m_row(row) {
                u = row.find_first();
            }

            Iterator(const AdjRow &row, Vertex start) : m_row(row), u(start) {}

            Vertex operator*() const { return u; }

            Iterator &operator++() {
                u = m_row.find_next(u);
                return *this;
            }

            bool operator==(const Iterator &other) const { return u == other.u; }

            bool operator!=(const Iterator &other) const { return !(*this == other); }
        };

        const AdjRow &m_row;
    public:
        explicit RowVertices(const AdjRow &row) : m_row(row) {}

        [[nodiscard]] Iterator begin() const { return Iterator(m_row); }

        [[nodiscard]] Iterator end() const { return Iterator(m_row, static_cast<Vertex>(AdjRow::npos)); }
    };

    static RowVertices vertices(const AdjRow &row) {
        return RowVertices(row);
    }

    [[nodiscard]] RowVertices neighbors(Vertex u) const {
        return RowVertices(adj[u]);
    }


    class VertexPairs {
        class Iterator {
            VertexPair m_uv;
            Vertex n;
        public:
            Iterator(VertexPair start, Vertex size) : m_uv(start), n(size) {}

            VertexPair operator*() const { return m_uv; }

            Iterator &operator++() {
                ++m_uv.v;
                if (m_uv.v == n) {
                    ++m_uv.u;
                    m_uv.v = m_uv.u + 1;
                }
                return *this;
            }

            bool operator==(const Iterator &other) const { return m_uv == other.m_uv; }

            bool operator!=(const Iterator &other) const { return !(*this == other); }
        };

        Vertex n;
    public:
        explicit VertexPairs(Vertex size) : n(size) {}

        [[nodiscard]] Iterator begin() const { return Iterator({0, 1}, n); }

        [[nodiscard]] Iterator end() const { return Iterator({n - 1, n}, n); }
    };

    [[nodiscard]] VertexPairs vertexPairs() const {
        return VertexPairs(n);
    }


    class Edges {
        class Iterator {
            const AdjMatrix &m_adj;
            VertexPair m_uv;
        public:
            Iterator(const AdjMatrix &adj, VertexPair start) : m_adj(adj), m_uv(start) {}

            explicit Iterator(const AdjMatrix &adj) : m_adj(adj), m_uv({0, 1}) {
                while (m_uv.u < adj.size() && m_adj[m_uv.u].none()) m_uv.u++;
                if (m_uv.u < adj.size())
                    m_uv.v = m_adj[m_uv.u].find_first();
            }

            VertexPair operator*() const { return m_uv; }

            Iterator &operator++() {
                m_uv.v = m_adj[m_uv.u].find_next(m_uv.v);
                while (m_uv.v >= m_adj.size() && m_uv.u < m_adj.size()) {
                    ++m_uv.u;
                    if (m_uv.u < m_adj.size())
                        m_uv.v = m_adj[m_uv.u].find_next(m_uv.u);
                }
                return *this;
            }

            bool operator==(const Iterator &other) const { return m_uv == other.m_uv; }

            bool operator!=(const Iterator &other) const { return !(*this == other); }

        };

        const AdjMatrix &m_adj;
    public:
        explicit Edges(const AdjMatrix &adj) : m_adj(adj) {}

        [[nodiscard]] Iterator begin() const { return Iterator(m_adj); }

        [[nodiscard]] Iterator end() const {
            return Iterator(m_adj, {static_cast<Vertex>(m_adj.size()), static_cast<Vertex>(AdjRow::npos)});
        }
    };

    [[nodiscard]] Edges edges() const {
        return Edges(adj);
    }


    /**
     * Iterate over all neighbors of u.
     * If the callback returns true, the iteration is stopped early.
     *
     * @tparam VertexCallback
     * @param u
     * @param callback
     * @return
     */
    template<typename VertexCallback>
    bool for_neighbors_of(Vertex u, VertexCallback callback) const {
        Vertex v = adj[u].find_first();
        while (v < n) {
            if (callback(v)) return true;
            v = adj[u].find_next(v);
        }
        return false;
    }

    /**
     * Iterate over all vertices u (0..n-1).
     * If the callback returns true, the iteration is stopped early.
     *
     * @tparam VertexCallback
     * @param callback
     * @return
     */
    template<typename VertexCallback>
    bool for_all_vertices(VertexCallback callback) const {
        for (Vertex u = 0; u < n; ++u) {
            if (callback(u)) return true;
        }
        return false;
    }

    /**
     * Iterate over all edges {u, v}.
     * If the callback returns true, the iteration is stopped early.
     *
     * @tparam VertexPairCallback
     * @param callback
     * @return
     */
    template<typename VertexPairCallback>
    bool for_all_edges(VertexPairCallback callback) const {
        return for_all_vertices([&](Vertex u) {
            Vertex v = adj[u].find_next(u);
            while (v < n) {
                if (callback(VertexPair(u, v))) return true;
                v = adj[u].find_next(v);
            }
            return false;
        });
    }

    /**
     * Iterate over all vertex pairs {u, v} ({0, 1}, {0, 2}, ..., {n-2, n-1}).
     * If the callback returns true, the iteration is stopped early.
     *
     * @tparam VertexPairCallback
     * @param callback
     * @return
     */
    template<typename VertexPairCallback>
    bool for_all_vertex_pairs(VertexPairCallback callback) const {
        for (Vertex u = 0; u < n; ++u) {
            for (Vertex v = u + 1; v < n; ++v) {
                if (callback(VertexPair(u, v))) return true;
            }
        }
        return false;
    }

    friend std::ostream &operator<<(std::ostream &os, const Graph &graph) {
        for (Vertex u = 0; u < graph.size(); ++u) {
            for (Vertex v = 0; v < graph.size(); ++v) {
                os << graph.adj[u][v] << " ";
            }
            os << "\n";
        }
        return os;
    }

private:

    [[nodiscard]] AdjRow all_vertices() const {
        return AdjRow(n, 1);
    }

    /**
     * Iterate over the vertex set indicated by row.
     *
     * @tparam VertexCallback
     * @param row
     * @param callback
     * @return
     */
    template<typename VertexCallback>
    static bool iterate(const AdjRow &row, VertexCallback callback) {
        Vertex u = row.find_first();
        while (u < row.size()) {
            if (callback(u)) return true;
            u = row.find_next(u);
        }
        return false;
    }

    template<typename VertexCallback>
    static bool iterate(const AdjRow &A, const AdjRow &B, VertexCallback callback) {
        return iterate(A, [&](Vertex a) {
            return iterate(B, [&](Vertex b) {
                return callback(a, b);
            });
        });
    }

    friend class FinderI;

    friend class Finder::NaiveP3;

    friend class Finder::CenterP3;

    friend class Finder::NaiveC4P4;

    friend class Finder::CenterC4P4;

    template<int length> friend
    class detail::Center;

    template<int length, bool with_cycles> friend
    class detail::CenterFinderImpl;
};


#endif //CONCEPT_GRAPH_H
