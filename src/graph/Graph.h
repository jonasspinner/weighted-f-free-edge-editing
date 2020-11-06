#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GRAPH_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GRAPH_H

#include <iostream>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <cmath>

#include "VertexPair.h"
#include "DynamicBitset.h"


class Graph {
public:
    using AdjRow = dynamic_bitset::DynamicBitset<>;
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

    Graph(const Graph &) = delete;

    Graph(Graph &&other) noexcept: m_size(other.m_size), m_adj(std::move(other.m_adj)) {}

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
        return Graph{size};
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
            using difference_type = std::make_signed_t<Vertex>;
            using pointer = void;
            using reference = Vertex;
            using iterator_category = std::random_access_iterator_tag;

            constexpr Iterator() noexcept {};

            constexpr explicit Iterator(Vertex start) noexcept: m_u(start) {}

            [[nodiscard]] constexpr value_type operator*() const noexcept {
                return m_u;
            }

            [[nodiscard]] constexpr value_type operator[](difference_type n) const noexcept {
                return m_u + static_cast<Vertex>(n);
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

            constexpr Iterator &operator+=(difference_type n) noexcept {
                m_u += static_cast<Vertex>(n);
                return *this;
            }

            constexpr Iterator &operator-=(difference_type n) noexcept {
                m_u -= static_cast<Vertex>(n);
                return *this;
            }

            [[nodiscard]] constexpr Iterator operator+(difference_type n) const noexcept {
                Iterator copy(*this);
                copy += n;
                return copy;
            }

            [[nodiscard]] friend constexpr Iterator operator+(difference_type n, const Iterator &it) noexcept {
                return it + n;
            }

            [[nodiscard]] constexpr Iterator operator-(difference_type n) const noexcept {
                Iterator copy(*this);
                copy -= n;
                return copy;
            }

            [[nodiscard]] constexpr difference_type operator-(const Iterator &other) const noexcept {
                return static_cast<difference_type>(m_u - other.m_u);
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
            using InnerIt = dynamic_bitset::IndexIterator<AdjRow::block_type>;
            InnerIt m_it{};
        public:
            using value_type = Vertex;
            using difference_type = std::ptrdiff_t;
            using pointer = void;
            using reference = Vertex;
            using iterator_category = std::forward_iterator_tag;

            constexpr Iterator() noexcept = default;

            explicit constexpr Iterator(const AdjRow &row) noexcept: m_it(row) {}

            struct end_tag {
            };

            constexpr Iterator(const AdjRow &row, end_tag) noexcept: m_it(row, row.size()) {}

            [[nodiscard]] constexpr value_type operator*() const noexcept {
                return static_cast<Vertex>(*m_it);
            }

            inline Iterator &operator++() noexcept {
                ++m_it;
                return *this;
            }

            [[nodiscard]] constexpr bool operator==(const Iterator &other) const noexcept {
                return m_it == other.m_it;
            }

            [[nodiscard]] constexpr bool operator!=(const Iterator &other) const noexcept {
                return !(*this == other);
            }

            [[nodiscard]] constexpr bool operator<(const Iterator &other) const noexcept {
                return m_it < other.m_it;
            }

            [[nodiscard]] constexpr bool operator<=(const Iterator &other) const noexcept {
                return m_it <= other.m_it;
            }

            [[nodiscard]] constexpr bool operator>(const Iterator &other) const noexcept {
                return m_it > other.m_it;
            }

            [[nodiscard]] constexpr bool operator>=(const Iterator &other) const noexcept {
                return m_it >= other.m_it;
            }
        };

    private:
        const AdjRow &m_row;

        friend class Graph;

    public:
        constexpr explicit RowVertices(const AdjRow &row) noexcept: m_row(row) {}

    public:
        using const_iterator = Iterator;

        [[nodiscard]] constexpr const_iterator begin() const noexcept {
            return Iterator{m_row};
        }

        [[nodiscard]] constexpr const_iterator end() const noexcept {
            return Iterator{m_row, Iterator::end_tag{}};
        }
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
            Vertex m_size{0};
        public:
            using value_type = VertexPair;
            using difference_type = std::ptrdiff_t;
            using pointer = void;
            using reference = const VertexPair &;
            using iterator_category = std::random_access_iterator_tag;

            constexpr Iterator() noexcept {};

            constexpr Iterator(VertexPair start, Vertex size) noexcept: m_uv(start), m_size(size) {}

            [[nodiscard]] constexpr value_type operator*() const noexcept { return m_uv; }

            [[nodiscard]] constexpr value_type operator[](difference_type n) const noexcept {
                return *(*this + n);
            }

            constexpr Iterator &operator++() noexcept {
                ++m_uv.v;
                if (m_uv.v == m_size) {
                    ++m_uv.u;
                    m_uv.v = m_uv.u + 1;
                }
                return *this;
            }

            [[nodiscard]] constexpr Iterator operator++(int) noexcept {
                Iterator copy(*this);
                ++*this;
                return copy;
            }

            constexpr Iterator &operator--() noexcept {
                --m_uv.v;
                if (m_uv.v == m_uv.u) {
                    --m_uv.u;
                    m_uv.v = m_size - 1;
                }
                return *this;
            }

            [[nodiscard]] constexpr Iterator operator--(int) noexcept {
                Iterator copy(*this);
                --*this;
                return copy;
            }

            constexpr Iterator &operator+=(difference_type n) noexcept {
                const auto new_v = m_uv.v + static_cast<Vertex>(n);

                if (m_uv.u < new_v && new_v < m_size) {
                    m_uv.v = new_v;
                } else {
                    auto index = pair_to_index(m_uv, m_size);
                    index += static_cast<unsigned long>(n);
                    m_uv = index_to_pair(index, m_size);
                }

                return *this;
            }

            constexpr Iterator &operator-=(difference_type n) noexcept {
                return *this += -n;
            }

            [[nodiscard]] constexpr Iterator operator+(difference_type n) const noexcept {
                Iterator copy(*this);
                copy += n;
                return copy;
            }

            [[nodiscard]] friend constexpr Iterator operator+(difference_type n, const Iterator &it) noexcept {
                return it + n;
            }

            [[nodiscard]] constexpr Iterator operator-(difference_type n) const noexcept {
                Iterator copy(*this);
                copy -= n;
                return copy;
            }

            [[nodiscard]] constexpr difference_type operator-(const Iterator &other) const noexcept {
                // This index is the reverse order. Therefore the operands are switched in the return statement.
                const auto idx = [&](const auto &uv) constexpr {
                    auto v_ = m_size - 1 - static_cast<difference_type>(uv.u);
                    auto u_ = m_size - 1 - static_cast<difference_type>(uv.v);
                    return v_ * (v_ - 1) / 2 + u_;
                };
                return idx(other.m_uv) - idx(m_uv);
            }

            [[nodiscard]] constexpr bool operator==(const Iterator &other) const noexcept {
                return m_uv == other.m_uv;
            }

            [[nodiscard]] constexpr bool operator!=(const Iterator &other) const noexcept {
                return !(*this == other);
            }

            [[nodiscard]] constexpr bool operator<(const Iterator &other) const noexcept {
                return m_uv < other.m_uv;
            }

            [[nodiscard]] constexpr bool operator<=(const Iterator &other) const noexcept {
                return m_uv <= other.m_uv;
            }

            [[nodiscard]] constexpr bool operator>(const Iterator &other) const noexcept {
                return m_uv > other.m_uv;
            }

            [[nodiscard]] constexpr bool operator>=(const Iterator &other) const noexcept {
                return m_uv >= other.m_uv;
            }

            /**
             * Reference: https://stackoverflow.com/a/27088560
             *
             * @param uv
             * @param size
             * @return
             */
            static constexpr unsigned long pair_to_index(VertexPair uv, unsigned int size) noexcept {
                const unsigned long i = uv.u, j = uv.v, n = size;
                const auto k = (n * (n - 1) / 2) - (n - i) * ((n - i) - 1) / 2 + j - i - 1;
                return k;
            }

            /**
             * Both pair_to_index and index_to_pair were experimentally verified for all vertex pairs for size up to
             * 10000.
             *
             * Reference: https://stackoverflow.com/a/27088560
             *
             * @param index
             * @param size
             * @return
             */
            static constexpr VertexPair index_to_pair(unsigned long index, unsigned int size) noexcept {
                const auto k = static_cast<long>(index);
                const auto n = static_cast<long>(size);
                const auto i = n - 2 - static_cast<long>(std::sqrt(-8 * k + 4 * n * (n - 1) - 7) / 2.0 - 0.5);
                const auto j = k + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2;
                return {static_cast<Vertex>(i), static_cast<Vertex>(j)};
            }
        };

    private:
        Vertex m_number_of_nodes;
    public:
        using const_iterator = Iterator;
        using size_type = std::size_t;

        constexpr explicit VertexPairs(Vertex size) noexcept: m_number_of_nodes(size) {}

        [[nodiscard]] constexpr const_iterator begin() const noexcept {
            return Iterator{{0, 1}, m_number_of_nodes};
        }

        [[nodiscard]] constexpr const_iterator end() const noexcept {
            const auto n = m_number_of_nodes;
            const auto one_after_last = n == 0 ? VertexPair{0, 1} : VertexPair{n - 1, n};
            return Iterator{one_after_last, m_number_of_nodes};
        }

        [[nodiscard]] constexpr size_type size() const noexcept {
            const auto n = static_cast<size_type>(m_number_of_nodes);
            return n * (n - 1) / 2;
        }

        [[nodiscard]] constexpr bool empty() const noexcept {
            return size() == 0;
        }
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
            static constexpr auto end_vertex = static_cast<Vertex>(AdjRow::npos);
            static_assert(end_vertex == std::numeric_limits<Vertex>::max());

            const AdjMatrix *m_adj{nullptr};
            VertexPair m_uv{0, end_vertex};  // default is end iterator for size <= 1
        public:
            using value_type = VertexPair;
            using difference_type = std::ptrdiff_t;
            using pointer = void;
            using reference = VertexPair;
            using iterator_category = std::forward_iterator_tag;

            constexpr Iterator() noexcept {};

            struct end_tag {
            };

            Iterator(const AdjMatrix &adj, end_tag) noexcept: m_adj(std::addressof(adj)) {
                const auto size = static_cast<Vertex>(m_adj->size());
                if (size >= 1) {
                    // The last possible edge is (n-2, n-1). Therefore, n-1 is the end position of u.
                    m_uv.u = size - 1;
                }
            }

            explicit Iterator(const AdjMatrix &adj) noexcept: m_adj(std::addressof(adj)) {
                const auto size = m_adj->size();
                if (size < 2)  // For tiny (0 or 1) graphs, this is already equal to the end iterator.
                    return;
                assert(size >= 2);

                m_uv.v = static_cast<Vertex>((*m_adj)[m_uv.u].find_first());
                while (m_uv.u <= size - 2 && m_uv.v == end_vertex) {
                    ++m_uv.u;
                    m_uv.v = static_cast<Vertex>((*m_adj)[m_uv.u].find_next(m_uv.u));
                }

                assert((m_uv.u <= size - 2 && m_uv.v <= size - 1) ||
                       (m_uv.u == size - 1 && m_uv.v == end_vertex));
            }

            [[nodiscard]] value_type operator*() const noexcept {
                assert(m_adj->size() >= 2);
                assert(m_uv.u <= m_adj->size() - 2 && m_uv.v <= m_adj->size() - 1);
                return m_uv;
            }

            Iterator &operator++() noexcept {
                const auto size = m_adj->size();
                assert(size >= 2);
                assert(m_uv.u <= size - 2 && m_uv.v <= size - 1);

                m_uv.v = static_cast<Vertex>((*m_adj)[m_uv.u].find_next(m_uv.v));
                while (m_uv.u <= size - 2 && m_uv.v == end_vertex) {
                    ++m_uv.u;
                    m_uv.v = static_cast<Vertex>((*m_adj)[m_uv.u].find_next(m_uv.u));
                }

                assert((m_uv.u <= size - 2 && m_uv.v <= size - 1) ||
                       (m_uv.u == size - 1 && m_uv.v == end_vertex));
                return *this;
            }

            [[nodiscard]] constexpr bool operator==(const Iterator &other) const noexcept {
                return m_uv == other.m_uv;
            }

            [[nodiscard]] constexpr bool operator!=(const Iterator &other) const noexcept {
                return !(*this == other);
            }

            [[nodiscard]] constexpr bool operator<(const Iterator &other) const noexcept {
                return m_uv < other.m_uv;
            }

            [[nodiscard]] constexpr bool operator<=(const Iterator &other) const noexcept {
                return m_uv <= other.m_uv;
            }

            [[nodiscard]] constexpr bool operator>(const Iterator &other) const noexcept {
                return m_uv > other.m_uv;
            }

            [[nodiscard]] constexpr bool operator>=(const Iterator &other) const noexcept {
                return m_uv >= other.m_uv;
            }

        };

        using const_iterator = Iterator;
        using size_type = std::size_t;

    private:
        const_iterator m_begin;
        const_iterator m_end;
    public:
        explicit Edges(const AdjMatrix &adj) noexcept:
                m_begin(adj),
                m_end(adj, Iterator::end_tag{}) {}

        [[nodiscard]] constexpr const_iterator begin() const noexcept { return m_begin; }

        [[nodiscard]] constexpr const_iterator end() const noexcept { return m_end; }

        [[nodiscard]] constexpr bool empty() const noexcept { return begin() == end(); }
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
        std::size_t num_edges = 0;
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
        AdjRow row(size());
        row.set();
        return row;
    }

    [[nodiscard]] const auto &adj(Vertex u) const {
        return m_adj[u];
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_GRAPH_H
