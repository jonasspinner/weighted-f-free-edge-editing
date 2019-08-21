//
// Created by jonas on 27.06.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GRAPH_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GRAPH_H

#include <iostream>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <yaml-cpp/yaml.h>

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
    class FindImpl;

    template<int length, bool with_cycles>
    class FindNearImpl;
}


class VertexPair {
public:
    Vertex u;
    Vertex v;

    VertexPair(Vertex x, Vertex y) : u(x < y ? x : y), v(x < y ? y : x) { assert(x != y); }

    friend std::ostream &operator<<(std::ostream &os, VertexPair uv) {
        return os << "{" << uv.u << ", " << uv.v << "}";
    }

    friend YAML::Emitter &operator<<(YAML::Emitter &out, VertexPair uv) {
        out << YAML::Flow;
        out << YAML::BeginSeq << uv.u << uv.v << YAML::EndSeq;
        return out;
    }

    friend bool operator==(const VertexPair &uv, const VertexPair &xy) {
        return (uv.u == xy.u) && (uv.v == xy.v);
    }

    friend bool operator<(const VertexPair &uv, const VertexPair &xy) {
        return uv.u < xy.u || (uv.u == xy.u && uv.v < xy.v);
    }
};


class Graph {
public:
    using AdjRow = boost::dynamic_bitset<>;
    using AdjMatrix = std::vector<AdjRow>;

private:
    unsigned int m_size;
    AdjMatrix m_adj;

public:

    explicit Graph(unsigned int size) : m_size(size), m_adj(m_size, AdjRow(m_size)) {}

    /**
     * Returns the number of vertices.
     *
     * @return
     */
    [[nodiscard]] Vertex size() const { return m_size; }


    void clear_edges() {
        for (auto &row: m_adj) { row.reset(); }
    }

    /**
     * Toggles the edge.
     *
     * @param edge
     */
    void toggle_edge(VertexPair edge) {
        const auto[u, v] = edge;
        m_adj[u].flip(v);
        m_adj[v].flip(u);
    }

    /**
     * Returns the degree of vertex u.
     *
     * @param u
     * @return
     */
    [[nodiscard]] size_t degree(Vertex u) const { return m_adj[u].count(); };

    /**
     * Checks whether the edge is in the graph.
     *
     * @param edge
     * @return
     */
    [[nodiscard]] bool has_edge(VertexPair edge) const { return m_adj[edge.u][edge.v]; }

    /**
     * Inserts the edge into the Graph.
     *
     * @param edge
     */
    void set_edge(VertexPair edge) {
        const auto[u, v] = edge;
        m_adj[u].set(v);
        m_adj[v].set(u);
    }

    /**
     * Inserts the edges into the Graph.
     *
     * @param edges
     */
    void set_edges(const std::vector<VertexPair> &edges) {
        for (VertexPair uv : edges) {
            assert(!has_edge(uv));
            set_edge(uv);
        }
    }

    /**
     * Removes the edge from the Graph.
     *
     * @param edge
     */
    void clear_edge(VertexPair edge) {
        const auto[u, v] = edge;
        m_adj[u].reset(v);
        m_adj[v].reset(u);
    }


    class Vertices {
        class Iterator {
            Vertex v;
        public:
            using value_type = Vertex;
            using iterator_category = std::forward_iterator_tag;

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
        return Vertices(m_size);
    }


    class RowVertices {
        class Iterator {
            const AdjRow &m_row;
            Vertex u;
        public:
            static const Vertex end_vertex = static_cast<Vertex>(AdjRow::npos);

            using value_type = Vertex;
            using difference_type = std::ptrdiff_t;
            using pointer = const Vertex*;
            using reference = const Vertex&;
            using iterator_category = std::forward_iterator_tag;

            explicit Iterator(const AdjRow &row) : m_row(row) {
                u = row.find_first();
                if (u >= m_row.size()) u = m_row.size();
            }

            Iterator(const AdjRow &row, Vertex start) : m_row(row), u(start) {}

            Vertex operator*() const {
                assert(u != end_vertex);
                assert(u < m_row.size());
                return u;
            }

            Iterator &operator++() {
                assert(u != end_vertex);
                assert(u < m_row.size());
                u = m_row.find_next(u);
                if (u >= m_row.size()) u = m_row.size();
                return *this;
            }

            bool operator==(const Iterator &other) const { return &m_row == &other.m_row && u == other.u; }

            bool operator!=(const Iterator &other) const { return !(*this == other); }
        };

        const AdjRow &m_row;
    public:
        explicit RowVertices(const AdjRow &row) : m_row(row) {}

        [[nodiscard]] Iterator begin() const { return Iterator(m_row); }

        [[nodiscard]] Iterator end() const { return Iterator(m_row, m_row.size()); }
    };

    static RowVertices iterate(const AdjRow &row) {
        return RowVertices(row);
    }

    [[nodiscard]] RowVertices neighbors(Vertex u) const {
        return RowVertices(m_adj[u]);
    }


    class VertexPairs {
        class Iterator {
            VertexPair m_uv;
            Vertex n;
        public:
            using value_type = VertexPair;
            using difference_type = std::ptrdiff_t;
            using pointer = const VertexPair*;
            using reference = const VertexPair&;
            using iterator_category = std::forward_iterator_tag;

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

            bool operator==(const Iterator &other) const { return m_uv == other.m_uv && n == other.n; }

            bool operator!=(const Iterator &other) const { return !(*this == other); }
        };

        Vertex n;
    public:
        explicit VertexPairs(Vertex size) : n(size) {}

        [[nodiscard]] Iterator begin() const { return Iterator({0, 1}, n); }

        [[nodiscard]] Iterator end() const {
            if (n > 0) {
                return Iterator({n - 1, n}, n);
            } else {
                return begin();
            }
        }
    };

    [[nodiscard]] VertexPairs vertexPairs() const {
        return VertexPairs(m_size);
    }


    class Edges {
        class Iterator {
            const AdjMatrix &m_adj;
            VertexPair m_uv;
        public:
            using value_type = VertexPair;
            using difference_type = std::ptrdiff_t;
            using pointer = const VertexPair*;
            using reference = const VertexPair&;
            using iterator_category = std::forward_iterator_tag;

            Iterator(const AdjMatrix &adj, VertexPair start) : m_adj(adj), m_uv(start) {}

            explicit Iterator(const AdjMatrix &adj) : m_adj(adj), m_uv({0, 1}) {
                const Vertex size = m_adj.size();
                while (m_uv.u < size && m_adj[m_uv.u].none()) m_uv.u++;
                if (m_uv.u < size) {
                    m_uv.v = m_adj[m_uv.u].find_first();
                    assert(m_uv.u < size - 1);
                    assert(m_uv.v < size);
                } else {
                    m_uv = {size - 1, size};
                }
            }

            VertexPair operator*() const {
                assert(m_uv.u < m_adj.size() - 1);
                assert(m_uv.v < m_adj.size());
                return m_uv;
            }

            Iterator &operator++() {
                const Vertex size = m_adj.size();
                assert(m_uv.u < size - 1);
                assert(m_uv.v < size);
                m_uv.v = m_adj[m_uv.u].find_next(m_uv.v);
                while (m_uv.u < size - 1 && m_uv.v >= size) {
                    ++m_uv.u;
                    m_uv.v = m_adj[m_uv.u].find_next(m_uv.u);
                }
                if (m_uv.u >= size - 1 || m_uv.v >= size){
                    m_uv = {size - 1, size};
                }
                return *this;
            }

            bool operator==(const Iterator &other) const { return &m_adj == &other.m_adj && m_uv == other.m_uv; }

            bool operator!=(const Iterator &other) const { return !(*this == other); }

        };

        const AdjMatrix &m_adj;
    public:
        explicit Edges(const AdjMatrix &adj) : m_adj(adj) {}

        [[nodiscard]] Iterator begin() const { return Iterator(m_adj); }

        [[nodiscard]] Iterator end() const {
            Vertex size = m_adj.size();
            return Iterator(m_adj, {size - 1, size});
        }
    };

    [[nodiscard]] Edges edges() const {
        return Edges(m_adj);
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
        Vertex v = m_adj[u].find_first();
        while (v < m_size) {
            if (callback(v)) return true;
            v = m_adj[u].find_next(v);
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
        for (Vertex u = 0; u < m_size; ++u) {
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
            Vertex v = m_adj[u].find_next(u);
            while (v < m_size) {
                if (callback(VertexPair(u, v))) return true;
                v = m_adj[u].find_next(v);
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
        for (Vertex u = 0; u < m_size; ++u) {
            for (Vertex v = u + 1; v < m_size; ++v) {
                if (callback(VertexPair(u, v))) return true;
            }
        }
        return false;
    }

    friend std::ostream &operator<<(std::ostream &os, const Graph &graph) {
        for (Vertex u = 0; u < graph.size(); ++u) {
            for (Vertex v = 0; v < graph.size(); ++v) {
                os << graph.m_adj[u][v] << " ";
            }
            os << "\n";
        }
        return os;
    }

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const Graph &graph) {
        unsigned num_edges = 0;
        for (Vertex u : graph.vertices())
            num_edges += graph.degree(u);
        num_edges /= 2;

        out << YAML::BeginMap;
        out << YAML::Key << "num_vertices" << YAML::Value << graph.size();
        out << YAML::Key << "num_edges" << YAML::Value << num_edges;
        out << YAML::Key << "adj";
        out << YAML::Value << YAML::BeginMap;
        for (Vertex u : graph.vertices()) {
            out << YAML::Key << u;
            out << YAML::Value << YAML::Flow << YAML::BeginSeq;
            for (Vertex v : graph.neighbors(u))
                out << v;
            out << YAML::EndSeq;
        }
        out << YAML::EndMap << YAML::EndMap;
        return out;
    }

private:

    [[nodiscard]] AdjRow all_vertices() const {
        return AdjRow(m_size, 1);
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

    template<int length, bool with_cycles>
    friend class detail::FindImpl;

    template<int length, bool with_cycles>
    friend class detail::FindNearImpl;
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_GRAPH_H
