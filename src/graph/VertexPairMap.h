#ifndef WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIRMAP_H
#define WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIRMAP_H


#include <iomanip>

#include "Graph.h"


template<typename T>
class VertexPairMap {
    using reference = typename std::vector<T>::reference;
    using const_reference = typename std::vector<T>::const_reference;

    Vertex m_size;
    std::vector<T> m_values;
public:

    /**
     * A map from vertex pairs to values T. {V \choose 2} \mapsto T.
     *
     * @param size The number of vertices. Vertices are indexed by 0..size-1.
     * @param initial An initial value. Default constructed if not given.
     */
    explicit VertexPairMap(Vertex size, T initial = T())
            : m_size(size),
              m_values(idx({m_size - 2, m_size - 1}) + 1 /*last valid index plus 1 is the size*/, initial) {}

    [[nodiscard]] const_reference operator[](VertexPair edge) const {
        assert(edge.u < m_size && edge.v < m_size);
        return m_values[idx(edge)];
    }

    [[nodiscard]] reference operator[](VertexPair edge) {
        assert(edge.u < m_size && edge.v < m_size);
        return m_values[idx(edge)];
    }

    [[nodiscard]] auto keys() const {
        return Graph::VertexPairs(m_size);
    }

    friend std::ostream &operator<<(std::ostream &os, const VertexPairMap &map) {
        for (Vertex u = 0; u < map.m_size; ++u) {
            os << "[ ";
            for (Vertex v = 0; v < map.m_size; ++v) {
                T value = (u != v) ? map[{u, v}] : T();
                os << std::fixed << std::setw(6) << std::setprecision(2) << value << " ";
            }
            os << "]\n";
        }
        return os;
    }


    friend YAML::Emitter &operator<<(YAML::Emitter &out, const VertexPairMap &map) {
        using namespace YAML;
        out << BeginMap;
        out << Key << "size";
        out << Value << map.m_size;
        out << Key << "values";
        out << Value << BeginSeq;
        for (Vertex u = 0; u < map.m_size; ++u) {
            out << Flow << BeginSeq;
            for (Vertex v = u + 1; v < map.m_size; ++v) {
                out << map[{u, v}];
            }
            out << EndSeq;
        }
        out << EndSeq << EndMap;
        return out;
    }

    [[nodiscard]] Vertex size() const { return m_size; }

private:
    /**
     * Returns the index for the given pair of vertices.
     *
     * For u, v \in 0..n-1 returns indices from 0 to n * (n - 1) / 2 - 1. The formula assumes that u < v.
     *
     * @param uv
     * @return
     */
    static inline size_t idx(VertexPair uv) {
        return uv.v * (uv.v - 1) / 2 + uv.u;
    }
};


template<>
class VertexPairMap<bool> {
    using AdjRow = boost::dynamic_bitset<>;
    using AdjMatrix = std::vector<AdjRow>;

    unsigned int m_size;
    AdjMatrix m_values;

public:
    using const_reference = AdjRow::const_reference;

    class reference {
    private:
        VertexPairMap<bool> &m_map;
        VertexPair m_uv;

        friend class VertexPairMap<bool>;

        reference(VertexPairMap<bool> &map, VertexPair uv) : m_map(map), m_uv(uv) {}
    public:
        reference &operator=(bool value) {
            m_map.m_values[m_uv.u][m_uv.v] = value;
            m_map.m_values[m_uv.v][m_uv.u] = value;
            return *this;
        }

        operator bool() const {
            return m_map.m_values[m_uv.u][m_uv.v];
        }
    };


    explicit VertexPairMap(Vertex size, bool initial = false) : m_size(size), m_values(m_size, AdjRow(m_size)) {
        if (initial) {
            for (auto &row : m_values) {
                row.set();
            }
        }
    }

    [[nodiscard]] const_reference operator[](VertexPair uv) const {
        assert(uv.u < m_size && uv.v < m_size);
        return m_values[uv.u][uv.v];
    }

    [[nodiscard]] reference operator[](VertexPair uv) {
        assert(uv.u < m_size && uv.v < m_size);
        return reference(*this, uv);
    }

    [[nodiscard]] auto keys() const {
        return Graph::VertexPairs(m_size);
    }

    friend std::ostream &operator<<(std::ostream &os, const VertexPairMap &map) {
        for (Vertex u = 0; u < map.m_size; ++u) {
            os << "[ ";
            for (Vertex v = 0; v < map.m_size; ++v) {
                bool value = u == v ? false : map[{u, v}];
                os << value << " ";
            }
            os << "]\n";
        }
        return os;
    }

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const VertexPairMap &map) {
        using namespace YAML;
        out << BeginMap;
        out << Key << "size";
        out << Value << map.m_size;
        out << Key << "values";
        out << Value << BeginSeq;
        for (Vertex u = 0; u < map.m_size; ++u) {
            out << Flow << BeginSeq;
            for (Vertex v = u + 1; v < map.m_size; ++v) {
                out << map[{u, v}];
            }
            out << EndSeq;
        }
        out << EndSeq << EndMap;
        return out;
    }

    [[nodiscard]] Vertex size() const { return m_size; }

    /**
     * Sets all entries to false.
     */
    void clear() {
        for (auto &row : m_values) {
            row.reset();
        }
    }
};

#endif //WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIRMAP_H
