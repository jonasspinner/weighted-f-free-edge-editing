//
// Created by jonas on 27.06.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GRAPH_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GRAPH_H

#include <iostream>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <yaml-cpp/yaml.h>

#include "VertexPair.h"


class Graph {
public:
    using AdjRow = boost::dynamic_bitset<>;
    using AdjMatrix = std::vector<AdjRow>;

private:
    unsigned int m_size;
    AdjMatrix m_adj;

public:

    /**
     * A undirected graph. Vertices are represented as numbers from 0 to size - 1.
     *
     * @param size
     */
    explicit Graph(unsigned int size) : m_size(size), m_adj(m_size, AdjRow(m_size)) {}

    Graph(const Graph&) = delete;

    Graph(Graph&& other) : m_size(other.m_size), m_adj(std::move(other.m_adj)) {}

    friend void swap(Graph &a, Graph &b) {
        using std::swap;
        swap(a.m_size, b.m_size);
        swap(a.m_adj, b.m_adj);
    }

    Graph copy() const {
        Graph other(m_size);
        other.m_adj = m_adj;
        return other;
    }

    static Graph from_vertex_pairs(unsigned int size, const std::vector<VertexPair> &vertexPairs ) {
        Graph new_graph(size);
        for (auto uv : vertexPairs) {
            new_graph.setEdge(uv);
        }
        return new_graph;
    }

    static Graph make_path(unsigned int size) {
        Graph path(size);
        for (Vertex u = 0; u + 1 < size; ++u) {
            path.setEdge({u, u + 1});
        }
        return path;
    }

    static Graph make_cycle(unsigned int size) {
        Graph cycle(size);
        for (Vertex u = 0; u + 1 < size; ++u) {
            cycle.setEdge({u, u + 1});
        }
        cycle.setEdge({0, size - 1});
        return cycle;
    }

    static Graph make_empty(unsigned int size) {
        return Graph(size);
    }

    /**
     * Returns the number of vertices.
     *
     * @return
     */
    [[nodiscard]] Vertex size() const { return m_size; }


    /**
     * Clears all edges.
     */
    void clearEdges() {
        for (auto &row: m_adj) { row.reset(); }
    }

    /**
     * Toggles the edge.
     *
     * @param edge
     */
    void toggleEdge(VertexPair edge) {
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
    [[nodiscard]] bool hasEdge(VertexPair edge) const { return m_adj[edge.u][edge.v]; }

    /**
     * Inserts the edge into the graph.
     *
     * @param edge
     */
    void setEdge(VertexPair edge) {
        const auto[u, v] = edge;
        m_adj[u].set(v);
        m_adj[v].set(u);
    }

    /**
     * Inserts the edges into the graph.
     *
     * @param edges
     */
    void setEdges(const std::vector<VertexPair> &edges) {
        for (VertexPair uv : edges) {
            assert(!hasEdge(uv));
            setEdge(uv);
        }
    }

    /**
     * Removes the edge from the graph.
     *
     * @param edge
     */
    void clearEdge(VertexPair edge) {
        const auto[u, v] = edge;
        m_adj[u].reset(v);
        m_adj[v].reset(u);
    }


    class Vertices {
    public:
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

    private:
        Vertex n;
    public:
        explicit Vertices(Vertex size) : n(size) {}

        [[nodiscard]] Iterator begin() const { return Iterator(0); }

        [[nodiscard]] Iterator end() const { return Iterator(n); }
    };

    /**
     * Returns a range for the vertices of the graph (from 0 to size - 1).
     *
     * @return
     */
    [[nodiscard]] Vertices vertices() const {
        return Vertices(m_size);
    }


    class RowVertices {
    public:
        class Iterator {
            const AdjRow *m_row;
            Vertex u;
        public:
#ifndef NDEBUG
            static const Vertex end_vertex = static_cast<Vertex>(AdjRow::npos);
#endif

            using value_type = Vertex;
            using difference_type = std::ptrdiff_t;
            using pointer = const Vertex *;
            using reference = const Vertex &;
            using iterator_category = std::forward_iterator_tag;

            explicit Iterator(const AdjRow *row) : m_row(row) {
                assert(m_row != nullptr);
                u = m_row->find_first();
                if (u >= m_row->size()) u = m_row->size();
            }

            Iterator(const AdjRow *row, Vertex start) : m_row(row), u(start) {
                assert(m_row != nullptr);
            }

            Vertex operator*() const {
                assert(u != end_vertex);
                assert(u < m_row->size());
                return u;
            }

            Iterator &operator++() {
                assert(u != end_vertex);
                assert(u < m_row->size());
                u = m_row->find_next(u);
                if (u >= m_row->size()) u = m_row->size();
                return *this;
            }

            bool operator==(const Iterator &other) const { return m_row == other.m_row && u == other.u; }

            bool operator!=(const Iterator &other) const { return !(*this == other); }
        };

    private:
        const AdjRow *m_row;
    public:
        explicit RowVertices(const AdjRow &row) : m_row(std::addressof(row)) {}

        [[nodiscard]] Iterator begin() const { return Iterator(m_row); }

        [[nodiscard]] Iterator end() const { return Iterator(m_row, m_row->size()); }
    };

    /**
     * Returns a range for the vertices indicated by the given adjacency row.
     * The vertices are ordered by their value.
     *
     * @param row
     * @return
     */
    [[nodiscard]] static RowVertices iterate(const AdjRow &row) {
        return RowVertices(row);
    }

    /**
     * Returns a range for the neighbors of the vertex u.
     * The vertices are ordered by their value.
     *
     * @param u
     * @return
     */
    [[nodiscard]] RowVertices neighbors(Vertex u) const {
        return RowVertices(m_adj[u]);
    }


    class VertexPairs {
    public:
        class Iterator {
            VertexPair m_uv;
            Vertex n;
        public:
            using value_type = VertexPair;
            using difference_type = std::ptrdiff_t;
            using pointer = const VertexPair *;
            using reference = const VertexPair &;
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

    private:
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

    /**
     * Returns a range for all unordered pairs of vertices of the graph.
     * The pairs are lexicographically ordered by their first and second vertex.
     *
     * @return
     */
    [[nodiscard]] VertexPairs vertexPairs() const {
        return VertexPairs(m_size);
    }


    class Edges {
    public:
        class Iterator {
            const AdjMatrix *m_adj;
            VertexPair m_uv;
        public:
            using value_type = VertexPair;
            using difference_type = std::ptrdiff_t;
            using pointer = const VertexPair *;
            using reference = const VertexPair &;
            using iterator_category = std::forward_iterator_tag;

            Iterator(const AdjMatrix *adj, VertexPair start) : m_adj(adj), m_uv(start) {
                assert(m_adj != nullptr);
            }

            explicit Iterator(const AdjMatrix *adj) : m_adj(adj), m_uv({0, 1}) {
                assert(m_adj != nullptr);
                const Vertex size = m_adj->size();
                while (m_uv.u < size && (*m_adj)[m_uv.u].none()) m_uv.u++;
                if (m_uv.u < size) {
                    m_uv.v = (*m_adj)[m_uv.u].find_first();
                    assert(m_uv.u < size - 1);
                    assert(m_uv.v < size);
                } else {
                    m_uv = {size - 1, size};
                }
            }

            VertexPair operator*() const {
                assert(m_uv.u < m_adj->size() - 1);
                assert(m_uv.v < m_adj->size());
                return m_uv;
            }

            Iterator &operator++() {
                const Vertex size = m_adj->size();
                assert(m_uv.u < size - 1);
                assert(m_uv.v < size);
                m_uv.v = (*m_adj)[m_uv.u].find_next(m_uv.v);
                while (m_uv.u < size - 1 && m_uv.v >= size) {
                    ++m_uv.u;
                    m_uv.v = (*m_adj)[m_uv.u].find_next(m_uv.u);
                }
                if (m_uv.u >= size - 1 || m_uv.v >= size) {
                    m_uv = {size - 1, size};
                }
                return *this;
            }

            bool operator==(const Iterator &other) const { return m_adj == other.m_adj && m_uv == other.m_uv; }

            bool operator!=(const Iterator &other) const { return !(*this == other); }

        };

    private:
        const AdjMatrix *m_adj;
    public:
        explicit Edges(const AdjMatrix &adj) : m_adj(std::addressof(adj)) {}

        [[nodiscard]] Iterator begin() const { return Iterator(m_adj); }

        [[nodiscard]] Iterator end() const {
            Vertex size = m_adj->size();
            return Iterator(m_adj, {size - 1, size});
        }
    };

    /**
     * Returns a range for the edges of the graph.
     * The pairs are lexicographically ordered by their first and second vertex.
     *
     * @return
     */
    [[nodiscard]] Edges edges() const {
        return Edges(m_adj);
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
        using namespace YAML;
        unsigned num_edges = 0;
        for (Vertex u : graph.vertices())
            num_edges += graph.degree(u);
        num_edges /= 2;

        out << BeginMap;
        out << Key << "num_vertices" << Value << graph.size();
        out << Key << "num_edges" << Value << num_edges;
        out << Key << "adj";
        out << Value << BeginMap;
        for (Vertex u : graph.vertices()) {
            out << Key << u;
            out << Value << Flow << BeginSeq;
            for (Vertex v : graph.neighbors(u))
                out << v;
            out << EndSeq;
        }
        out << EndMap << EndMap;
        return out;
    }

    /**
     * Returns an adjacency row with every vertex.
     *
     * @return
     */
    [[nodiscard]] AdjRow full_adjacency_row() const {
        auto row = AdjRow(m_size);
        row.set();
        return row;
    }

    friend class FinderI;
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_GRAPH_H
