#ifndef WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIRMAP_H
#define WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIRMAP_H


#include <iomanip>

#include <boost/dynamic_bitset.hpp>

#include "Graph.h"


template<typename T>
class VertexPairMap {
    Vertex m_size;
    std::vector<T> m_values;
public:
    using reference = typename std::vector<T>::reference;
    using const_reference = typename std::vector<T>::const_reference;
    using value_type = T;

    /**
     * A map from vertex pairs to values T. {V \choose 2} \mapsto T.
     *
     * @param size The number of vertices. Vertices are indexed by 0..size-1.
     * @param initial An initial value. Default constructed if not given.
     */
    explicit VertexPairMap(Vertex size, T initial = T()) noexcept:
            m_size(size),
            m_values(idx({m_size - 2, m_size - 1}) + 1 /*last valid index plus 1 is the size*/, initial) {}

    [[nodiscard]] constexpr const_reference operator[](VertexPair uv) const noexcept {
        assert(uv.u < m_size && uv.v < m_size);
        return m_values[idx(uv)];
    }

    [[nodiscard]] constexpr reference operator[](VertexPair uv) noexcept {
        assert(uv.u < m_size && uv.v < m_size);
        return m_values[idx(uv)];
    }

    [[nodiscard]] constexpr auto keys() const noexcept {
        return Graph::VertexPairs{m_size};
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

    [[nodiscard]] constexpr auto size() const noexcept { return m_size; }

private:
    /**
     * Returns the index for the given pair of vertices.
     *
     * For u, v \in 0..n-1 returns indices from 0 to n * (n - 1) / 2 - 1. The formula assumes that u < v.
     *
     * @param uv
     * @return
     */
    static constexpr std::size_t idx(VertexPair uv) noexcept {
        auto u = static_cast<std::size_t>(uv.u);
        auto v = static_cast<std::size_t>(uv.v);
        return v * (v - 1) / 2 + u;
    }
};


template<>
class VertexPairMap<bool> {
    using Values = boost::dynamic_bitset<>;

    unsigned int m_size;
    Values m_values;

public:
    using reference = Values::reference;
    using const_reference = Values::const_reference;
    using value_type = bool;

    explicit VertexPairMap(Vertex size, bool initial = false) noexcept:
            m_size(size),
            m_values(idx({m_size - 2, m_size - 1}) + 1) {
        if (initial)
            m_values.set();
    }

    [[nodiscard]] const_reference operator[](VertexPair uv) const noexcept {
        assert(uv.u < m_size && uv.v < m_size);
        return m_values[idx(uv)];
    }

    [[nodiscard]] reference operator[](VertexPair uv) noexcept {
        assert(uv.u < m_size && uv.v < m_size);
        return m_values[idx(uv)];
    }

    [[nodiscard]] constexpr auto keys() const noexcept {
        return Graph::VertexPairs{m_size};
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

    [[nodiscard]] constexpr auto size() const noexcept { return m_size; }

    /**
     * Sets all entries to false.
     */
    void clear() noexcept {
        m_values.reset();
    }

    /**
 * Returns the index for the given pair of vertices.
 *
 * For u, v \in 0..n-1 returns indices from 0 to n * (n - 1) / 2 - 1. The formula assumes that u < v.
 *
 * @param uv
 * @return
 */
    static constexpr std::size_t idx(VertexPair uv) noexcept {
        auto u = static_cast<std::size_t>(uv.u);
        auto v = static_cast<std::size_t>(uv.v);
        return v * (v - 1) / 2 + u;
    }
};

#endif //WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIRMAP_H
