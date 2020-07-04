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
    unsigned int m_size{0};
    AdjMatrix m_adj;

public:

    /**
     * A undirected graph. Vertices are consecutive integers from 0 to size - 1.
     *
     * @param size
     */
    explicit Graph(unsigned int size) noexcept: m_size(size), m_adj(m_size, AdjRow(m_size)) {}

    Graph(const Graph&) = delete;

    Graph(Graph&& other) noexcept: m_size(other.m_size), m_adj(std::move(other.m_adj)) {}

    friend void swap(Graph &a, Graph &b) noexcept {
        using std::swap;
        swap(a.m_size, b.m_size);
        swap(a.m_adj, b.m_adj);
    }

    Graph copy() const noexcept {
        Graph other(size());
        other.m_adj = m_adj;
        return other;
    }

    static Graph from_edges(unsigned int size, const std::vector<VertexPair> &edges) noexcept {
        Graph new_graph(size);
        for (auto uv : edges) {
            new_graph.setEdge(uv);
        }
        return new_graph;
    }

    static Graph make_path_graph(unsigned int size) noexcept {
        Graph path(size);
        for (Vertex u = 0; u + 1 < size; ++u) {
            path.setEdge({u, u + 1});
        }
        return path;
    }

    static Graph make_cycle_graph(unsigned int size) noexcept {
        Graph cycle(size);
        for (Vertex u = 0; u + 1 < size; ++u) {
            cycle.setEdge({u, u + 1});
        }
        cycle.setEdge({0, size - 1});
        return cycle;
    }

    static Graph make_empty_graph(unsigned int size) noexcept {
        return Graph(size);
    }

    /**
     * Returns the number of vertices.
     *
     * @return
     */
    [[nodiscard]] constexpr unsigned int size() const noexcept { return m_size; }

    [[nodiscard]] constexpr unsigned int number_of_vertices() const noexcept { return size(); }


    /**
     * Clears all edges.
     */
    void clearEdges() noexcept {
        for (auto &row: m_adj) {
            row.reset();
        }
    }

    /**
     * Toggles the edge.
     *
     * @param edge
     */
    inline void toggleEdge(VertexPair edge) noexcept {
        const auto[u, v] = edge;
        assert(u < size());
        assert(v < size());
        m_adj[u].flip(v);
        m_adj[v].flip(u);
    }

    /**
     * Returns the degree of vertex u.
     *
     * @param u
     * @return
     */
    [[nodiscard]] inline size_t degree(Vertex u) const noexcept {
        assert(u < size());
        return m_adj[u].count();
    };

    /**
     * Checks whether the edge is in the graph.
     *
     * @param edge
     * @return
     */
    [[nodiscard]] inline bool hasEdge(VertexPair edge) const noexcept {
        assert(edge.u < size());
        assert(edge.v < size());
        return m_adj[edge.u][edge.v];
    }

    /**
     * Inserts the edge into the graph.
     *
     * @param edge
     */
    void inline setEdge(VertexPair edge) noexcept {
        const auto[u, v] = edge;
        assert(u < size());
        assert(v < size());
        m_adj[u].set(v);
        m_adj[v].set(u);
    }

    /**
     * Removes the edge from the graph.
     *
     * @param edge
     */
    void clearEdge(VertexPair edge) noexcept {
        const auto[u, v] = edge;
        assert(u < size());
        assert(v < size());
        m_adj[u].reset(v);
        m_adj[v].reset(u);
    }


    class Vertices {
    public:
        class Iterator {
            Vertex m_u{0};
        public:
            using value_type = Vertex;
            using difference_type = Vertex;
            using pointer = void;
            using reference = Vertex;
            using iterator_category = std::random_access_iterator_tag;

            constexpr Iterator() noexcept {};

            constexpr explicit Iterator(Vertex start) noexcept: m_u(start) {}

            [[nodiscard]] constexpr value_type operator*() const noexcept {
                return m_u;
            }

            [[nodiscard]] constexpr value_type operator[](std::size_t n) const noexcept {
                return *(*this + n);
            }

            constexpr Iterator &operator++() noexcept {
                ++m_u;
                return *this;
            }

            [[nodiscard]] constexpr Iterator operator++(int) noexcept {
                Iterator copy(*this);
                ++m_u;
                return copy;
            }

            constexpr Iterator &operator--() noexcept {
                --m_u;
                return *this;
            }

            [[nodiscard]] constexpr Iterator operator--(int) noexcept {
                Iterator copy(*this);
                --m_u;
                return copy;
            }

            constexpr Iterator &operator+=(unsigned int n) noexcept {
                m_u += n;
                return *this;
            }

            constexpr Iterator &operator-=(unsigned int n) noexcept {
                m_u -= n;
                return *this;
            }

            [[nodiscard]] constexpr Iterator operator+(unsigned int n) const noexcept {
                Iterator copy(*this);
                copy += n;
                return copy;
            }

            [[nodiscard]] friend constexpr auto operator+(unsigned int n, const Iterator &it) noexcept {
                return it + n;
            }

            [[nodiscard]] constexpr Iterator operator-(unsigned int n) const noexcept {
                Iterator copy(*this);
                copy -= n;
                return copy;
            }

            [[nodiscard]] constexpr difference_type operator-(const Iterator &other) const noexcept {
                return m_u - other.m_u;
            }

            [[nodiscard]] constexpr bool operator==(const Iterator &other) const noexcept {
                return m_u == other.m_u;
            }

            [[nodiscard]] constexpr bool operator!=(const Iterator &other) const noexcept {
                return !(*this == other);
            }

            [[nodiscard]] constexpr bool operator<(const Iterator &other) const noexcept {
                return m_u < other.m_u;
            }

            [[nodiscard]] constexpr bool operator<=(const Iterator &other) const noexcept {
                return m_u <= other.m_u;
            }

            [[nodiscard]] constexpr bool operator>(const Iterator &other) const noexcept {
                return m_u > other.m_u;
            }

            [[nodiscard]] constexpr bool operator>=(const Iterator &other) const noexcept {
                return m_u >= other.m_u;
            }
        };

    private:
        Vertex m_number_of_nodes;
    public:
        using const_iterator = Iterator;
        using size_type = std::size_t;

        constexpr explicit Vertices(Vertex size) noexcept: m_number_of_nodes(size) {}

        [[nodiscard]] constexpr const_iterator begin() const noexcept { return Iterator{}; }

        [[nodiscard]] constexpr const_iterator end() const noexcept { return Iterator{m_number_of_nodes}; }

        [[nodiscard]] constexpr size_type size() const noexcept { return m_number_of_nodes; }

        [[nodiscard]] constexpr bool empty() const noexcept { return size() == 0; }
    };

    /**
     * Returns a range for the vertices of the graph (from 0 to size - 1).
     *
     * @return
     */
    [[nodiscard]] constexpr auto vertices() const noexcept {
        return Vertices{size()};
    }


    class RowVertices {
    public:
        class Iterator {
            const AdjRow *m_row{nullptr};
            Vertex m_u{0};
        public:
#ifndef NDEBUG
            static const Vertex end_vertex = static_cast<Vertex>(AdjRow::npos);
#endif

            using value_type = Vertex;
            using difference_type = std::ptrdiff_t;
            using pointer = void;
            using reference = Vertex;
            using iterator_category = std::forward_iterator_tag;

            constexpr Iterator() noexcept {};

            constexpr explicit Iterator(const AdjRow *row) noexcept: m_row(row) {
                assert(m_row != nullptr);
                m_u = m_row->find_first();
            }

            constexpr Iterator(const AdjRow *row, Vertex start) noexcept: m_row(row), m_u(start) {
                assert(m_row != nullptr);
            }

            Vertex operator*() const noexcept {
                assert(m_u != end_vertex);
                assert(m_u < m_row->size());
                return m_u;
            }

            Iterator &operator++() noexcept {
                assert(m_u != end_vertex);
                assert(m_u < m_row->size());
                m_u = m_row->find_next(m_u);
                return *this;
            }

            constexpr bool operator==(const Iterator &other) const noexcept {
                return m_row == other.m_row && m_u == other.m_u;
            }

            constexpr bool operator!=(const Iterator &other) const noexcept { return !(*this == other); }

            [[nodiscard]] constexpr bool operator<(const Iterator &other) const noexcept {
                return m_u < other.m_u;
            }

            [[nodiscard]] constexpr bool operator<=(const Iterator &other) const noexcept {
                return m_u <= other.m_u;
            }

            [[nodiscard]] constexpr bool operator>(const Iterator &other) const noexcept {
                return m_u > other.m_u;
            }

            [[nodiscard]] constexpr bool operator>=(const Iterator &other) const noexcept {
                return m_u >= other.m_u;
            }
        };

    private:
        const AdjRow *m_row;

        friend class Graph;
        constexpr explicit RowVertices(const AdjRow &row) noexcept: m_row(std::addressof(row)) {}
    public:
        using const_iterator = Iterator;

        [[nodiscard]] constexpr const_iterator begin() const noexcept { return Iterator{m_row}; }

        [[nodiscard]] constexpr const_iterator end() const noexcept {
            return Iterator{m_row, static_cast<Vertex>(AdjRow::npos)};
        }

        [[nodiscard]] constexpr bool empty() const noexcept { return begin() == end(); }
    };

    /**
     * Returns a range for the vertices indicated by the given adjacency row.
     * The vertices are ordered by their value.
     *
     * @param row
     * @return
     */
    [[nodiscard]] static constexpr auto iterate(const AdjRow &row) noexcept {
        return RowVertices{row};
    }

    /**
     * Returns a range for the neighbors of the vertex u.
     * The vertices are ordered by their value.
     *
     * @param u
     * @return
     */
    [[nodiscard]] auto neighbors(Vertex u) const noexcept {
        return RowVertices{m_adj[u]};
    }


    class VertexPairs {
    public:
        class Iterator {
            VertexPair m_uv{};
            Vertex n{0};
        public:
            using value_type = VertexPair;
            using difference_type = std::ptrdiff_t;
            using pointer = const VertexPair *;
            using reference = const VertexPair &;
            using iterator_category = std::forward_iterator_tag;

            constexpr Iterator() noexcept {};

            constexpr Iterator(VertexPair start, Vertex size) noexcept: m_uv(start), n(size) {}

            constexpr VertexPair operator*() const noexcept { return m_uv; }

            constexpr Iterator &operator++() noexcept {
                ++m_uv.v;
                if (m_uv.v == n) {
                    ++m_uv.u;
                    m_uv.v = m_uv.u + 1;
                }
                return *this;
            }

            [[nodiscard]] constexpr bool operator==(const Iterator &other) const noexcept {
                return m_uv == other.m_uv && n == other.n;
            }

            [[nodiscard]] constexpr bool operator!=(const Iterator &other) const noexcept { return !(*this == other); }
        };

    private:
        Vertex m_number_of_nodes;
    public:
        using const_iterator = Iterator;
        using size_type = std::size_t;

        constexpr explicit VertexPairs(Vertex size) noexcept: m_number_of_nodes(size) {}

        [[nodiscard]] constexpr Iterator begin() const noexcept {
            return Iterator({0, 1}, m_number_of_nodes);
        }

        [[nodiscard]] constexpr Iterator end() const noexcept {
            return m_number_of_nodes == 0 ?
                begin() :
                Iterator({m_number_of_nodes - 1, m_number_of_nodes}, m_number_of_nodes);
        }

        [[nodiscard]] constexpr size_type size() const noexcept {
            return m_number_of_nodes * (m_number_of_nodes - 1) / 2;
        }

        [[nodiscard]] constexpr bool empty() const noexcept { return size() == 0; }
    };

    /**
     * Returns a range for all unordered pairs of vertices of the graph.
     * The pairs are lexicographically ordered by their first and second vertex.
     *
     * @return
     */
    [[nodiscard]] constexpr auto vertexPairs() const noexcept {
        return VertexPairs{size()};
    }


    class Edges {
    public:
        class Iterator {
            const AdjMatrix *m_adj{nullptr};
            VertexPair m_uv{};
        public:
            using value_type = VertexPair;
            using difference_type = std::ptrdiff_t;
            using pointer = void;
            using reference = VertexPair;
            using iterator_category = std::forward_iterator_tag;

            constexpr Iterator() noexcept {};

            constexpr Iterator(const AdjMatrix *adj, VertexPair start) noexcept: m_adj(adj), m_uv(start) {
                assert(m_adj != nullptr);
            }

            explicit Iterator(const AdjMatrix *adj) noexcept: m_adj(adj), m_uv({0, 1}) {
                assert(m_adj != nullptr);
                const Vertex size = m_adj->size();
                while (m_uv.u < size && (*m_adj)[m_uv.u].none()) {
                    ++m_uv.u;
                }
                if (m_uv.u < size) {
                    m_uv.v = (*m_adj)[m_uv.u].find_first();
                    assert(m_uv.u < size - 1);
                    assert(m_uv.v < size);
                } else {
                    m_uv = {size - 1, size};
                }
            }

            value_type operator*() const noexcept {
                assert(m_adj->empty() || m_uv.u < m_adj->size() - 1);
                assert(m_uv.v < m_adj->size());
                return m_uv;
            }

            Iterator &operator++() noexcept {
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

            [[nodiscard]] constexpr bool operator==(const Iterator &other) const noexcept {
                return m_adj == other.m_adj && m_uv == other.m_uv;
            }

            [[nodiscard]] constexpr bool operator!=(const Iterator &other) const noexcept { return !(*this == other); }

        };

        using const_iterator = Iterator;
        using size_type = std::size_t;

    private:
        const_iterator m_begin;
        const_iterator m_end;
    public:
        explicit Edges(const AdjMatrix &adj) noexcept:
                m_begin(Iterator{std::addressof(adj)}),
                m_end(Iterator{std::addressof(adj), {static_cast<Vertex>(adj.size() - 1),
                                       static_cast<Vertex>(adj.size())}}) {}

        [[nodiscard]] const_iterator begin() const noexcept { return m_begin; }

        [[nodiscard]] const_iterator end() const noexcept { return m_end; }

        [[nodiscard]] bool empty() const noexcept { return begin() == end(); }
    };

    /**
     * Returns a range for the edges of the graph.
     * The pairs are lexicographically ordered by their first and second vertex.
     *
     * @return
     */
    [[nodiscard]] auto edges() const noexcept {
        return Edges{m_adj};
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
    [[nodiscard]] AdjRow full_adjacency_row() const noexcept {
        AdjRow row(m_size);
        row.set();
        return row;
    }

    friend class FinderI;
    friend class CenterC4P4Finder;
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_GRAPH_H
