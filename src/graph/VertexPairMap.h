//
// Created by jonas on 02.07.19.
//

#ifndef CONCEPT_VERTEXPAIRMAP_H
#define CONCEPT_VERTEXPAIRMAP_H


#include <iomanip>

#include "Graph.h"

template<typename T>
class VertexPairMap {
    using reference = typename std::vector<T>::reference;
    using const_reference = typename std::vector<T>::const_reference;

    Vertex n;
    std::vector<T> values;
public:

    explicit VertexPairMap(Vertex size, T initial = T()) : n(size), values(idx({n - 2, n - 1}) + 1, initial) {}

    const_reference operator[](VertexPair edge) const {
        assert(edge.u < n && edge.v < n);
        return values[idx(edge)];
    }

    reference operator[](VertexPair edge) {
        assert(edge.u < n && edge.v < n);
        return values[idx(edge)];
    }

    friend std::ostream &operator<<(std::ostream &os, const VertexPairMap &map) {
        for (Vertex u = 0; u < map.n; ++u) {
            os << "[ ";
            for (Vertex v = 0; v < map.n; ++v) {
                T value = (u != v) ? map[{u, v}] : T();
                os << std::fixed << std::setw(6) << std::setprecision(2) << value << " ";
            }
            os << "]\n";
        }
        return os;
    }


    friend YAML::Emitter &operator<<(YAML::Emitter &out, const VertexPairMap &map) {
        out << YAML::BeginMap;
        out << YAML::Key << "size";
        out << YAML::Value << map.n;
        out << YAML::Key << "values";
        out << YAML::Value << YAML::BeginSeq;
        for (Vertex u = 0; u < map.n; ++u) {
            out << YAML::Flow << YAML::BeginSeq;
            for (Vertex v = u + 1; v < map.n; ++v) {
                out << map[{u, v}];
            }
            out << YAML::EndSeq;
        }
        out << YAML::EndSeq << YAML::EndMap;
        return out;
    }

    [[nodiscard]] Vertex size() const { return n; }

private:
    static inline size_t idx(VertexPair uv) {
        return uv.v * (uv.v - 1) / 2 + uv.u;
    }
};


template<>
class VertexPairMap<bool> {
    using Block = uint_fast32_t;
    using AdjRow = boost::dynamic_bitset<Block>;
    using AdjMatrix = std::vector<AdjRow>;

    unsigned int n;
    AdjMatrix adj;

public:
    using const_reference = AdjRow::const_reference;

    class reference {
        VertexPairMap<bool> &m_map;
        VertexPair m_uv;
    public:
        reference(VertexPairMap<bool> &map, VertexPair uv) : m_map(map), m_uv(uv) {}

        reference &operator=(bool value) {
            m_map.adj[m_uv.u][m_uv.v] = value;
            m_map.adj[m_uv.v][m_uv.u] = value;
            return *this;
        }

        operator bool() const {
            return m_map.adj[m_uv.u][m_uv.v];
        }
    };


    explicit VertexPairMap(Vertex size, bool initial = false) : n(size),
                                                             adj(n, AdjRow(n)) {
        if (initial) {
            for (auto &row : adj) {
                row.set();
            }
        }
    }

    const_reference operator[](VertexPair uv) const {
        assert(uv.u < n && uv.v < n);
        return adj[uv.u][uv.v];
    }

    reference operator[](VertexPair uv) {
        assert(uv.u < n && uv.v < n);
        return reference(*this, uv);
    }

    [[nodiscard]] Vertex size() const { return n; }
};

#endif //CONCEPT_VERTEXPAIRMAP_H
