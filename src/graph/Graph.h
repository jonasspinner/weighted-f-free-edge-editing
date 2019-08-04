//
// Created by jonas on 27.06.19.
//

#ifndef CONCEPT_GRAPH_H
#define CONCEPT_GRAPH_H

#include <iostream>
#include <vector>
#include <boost/dynamic_bitset.hpp>

#include "../Definitions.h"

namespace Finder {
    class NaiveC4P4;
    class CenterC4P4;
    class CenterP3;
    class NaiveP3;
}
// template <int length> class Center;
namespace detail {
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
    /*
    using const_reference = boost::dynamic_bitset<Block>::const_reference;
    class reference {
        Graph& graph; VertexPair uv;
    public:
        reference(Graph& graph, VertexPair uv) : graph(graph), uv(uv) {}
        reference& operator=(bool value) {
            if (value) graph.set_edge(uv); else graph.clear_edge(uv);
            return *this;
        }
        explicit operator bool() const {
            return graph.adj[uv.u][uv.v];
        }
    };*/

    explicit Graph(unsigned int n) : n(n), adj(n, boost::dynamic_bitset<Block>(n)) {}

    /*
    const_reference operator[](VertexPair uv) const {
        return adj[uv.u][uv.v];
    }

    reference operator[](VertexPair uv) {
        return reference(*this, uv);
    } */

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
    template <int length, bool with_cycles> friend class detail::CenterFinderImpl;
    // template <int length> friend class Center;
};


#endif //CONCEPT_GRAPH_H
