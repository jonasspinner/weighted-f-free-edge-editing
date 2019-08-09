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
    template <int length> class Center;
    template <int length, bool with_cycles> class CenterFinderImpl;
}


class VertexPair {
public:
    Vertex u;
    Vertex v;

    VertexPair(Vertex x, Vertex y) : u(x < y ? x : y), v(x < y ? y : x) { assert(x != y); }

    friend std::ostream &operator<<(std::ostream &os, const VertexPair &uv) {
        return os << "{" << uv.u << ", " << uv.v << "}";
    }

    friend bool operator==(const VertexPair& uv, const VertexPair &xy) {
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

    explicit Graph(unsigned int n) : n(n), adj(n, boost::dynamic_bitset<Block>(n)) {}

    /**
     * Returns the number of vertices.
     *
     * @return
     */
    [[nodiscard]] Vertex size() const { return n; }


    void clear_edges() {
        for (auto& row: adj) { row.reset(); }
    }

    /**
     * Toggles the edge.
     *
     * @param edge
     */
    void toggle_edge(const VertexPair &edge) {
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
    [[nodiscard]] bool has_edge(const VertexPair &edge) const { return adj[edge.u][edge.v]; }

    /**
     * Sets the edge.
     *
     * @param edge
     */
    void set_edge(const VertexPair &edge) {
        const auto&[u, v] = edge;
        adj[u][v] = true;
        adj[v][u] = true;
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
    void clear_edge(const VertexPair &edge) {
        const auto&[u, v] = edge;
        adj[u][v] = false;
        adj[v][u] = false;
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
            bool operator==(const Iterator& other) { return v == other.v; }
            bool operator!=(const Iterator& other) { return !(*this == other); }
        };
        Vertex n;
    public:
        explicit Vertices(Vertex n) : n(n) {}
        [[nodiscard]] Iterator begin() const { return Iterator(0); }
        [[nodiscard]] Iterator end() const { return Iterator(n); }
    };
    [[nodiscard]] Vertices vertices() const {
        return Vertices(n);
    }


    class RowVertices {
        class Iterator {
            Vertex u;
            const AdjRow &row;
        public:
            explicit Iterator(const AdjRow &row) : row(row) {
                u = row.find_first();
            }
            Iterator(const AdjRow &row, Vertex start) : row(row), u(start) {}
            Vertex operator*() const { return u; }
            Iterator &operator++() {
                u = row.find_next(u);
                return *this;
            }
            bool operator==(const Iterator& other) { return u == other.u; }
            bool operator!=(const Iterator& other) { return !(*this == other); }
        };
        const AdjRow &row;
    public:
        explicit RowVertices(const AdjRow &row) : row(row) {}
        [[nodiscard]] Iterator begin() const { return Iterator(row); }
        [[nodiscard]] Iterator end() const { return Iterator(row, static_cast<Vertex>(AdjRow::npos)); }
    };
    static RowVertices vertices(const AdjRow &row) {
        return RowVertices(row);
    }

    [[nodiscard]] RowVertices neighbors(Vertex u) const {
        return RowVertices(adj[u]);
    }


    class VertexPairs {
        class Iterator {
            VertexPair uv;
            Vertex n;
        public:
            Iterator(VertexPair start, Vertex n) : uv(start), n(n) { }
            VertexPair operator*() const { return uv; }
            Iterator &operator++() {
                ++uv.v;
                if (uv.v == n) { ++uv.u; uv.v = uv.u + 1; }
                return *this;
            }
            bool operator==(const Iterator& other) { return uv == other.uv; }
            bool operator!=(const Iterator& other) { return !(*this == other); }
        };
        Vertex n;
    public:
        explicit VertexPairs(Vertex n) : n(n) {}
        [[nodiscard]] Iterator begin() const { return Iterator({0, 1}, n); }
        [[nodiscard]] Iterator end() const { return Iterator({n - 1, n}, n); }
    };
    [[nodiscard]] VertexPairs vertexPairs() const {
        return VertexPairs(n);
    }


    class Edges {
        class Iterator {
            const AdjMatrix &adj;
            VertexPair uv;
        public:
            Iterator(const AdjMatrix &adj, VertexPair start) : adj(adj), uv(start) {}
            explicit Iterator(const AdjMatrix &adj) : adj(adj), uv({0, 1}) {
                while (adj[uv.u].count() == 0) uv.u++;
                uv.v = adj[uv.u].find_first();
            }
            VertexPair operator*() const { return uv; }
            Iterator &operator++() {
                uv.v = adj[uv.u].find_next(uv.v);
                while (uv.v >= adj.size() && uv.u < adj.size()) {
                    ++uv.u;
                    if (uv.u < adj.size())
                        uv.v = adj[uv.u].find_next(uv.u);
                }
                return *this;
            }
            bool operator==(const Iterator& other) { return uv == other.uv; }
            bool operator!=(const Iterator& other) { return !(*this == other); }

        };
        const AdjMatrix& adj;
    public:
        explicit Edges(const AdjMatrix& adj) : adj(adj) {}
        [[nodiscard]] Iterator begin() const { return Iterator(adj); }
        [[nodiscard]] Iterator end() const { return Iterator(adj, {static_cast<Vertex>(adj.size()), static_cast<Vertex>(AdjRow::npos)}); }
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
        return for_all_vertices([&](auto u) {
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

    template <typename VertexCallback>
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
    template <int length> friend class detail::Center;
    template <int length, bool with_cycles> friend class detail::CenterFinderImpl;
};


#endif //CONCEPT_GRAPH_H
