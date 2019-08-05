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
    // explicit VertexPairMap(const Graph &graph, T initial = T()) : n(graph.size()),
    //                                                               values(idx({n - 2, n - 1}) + 1, initial) {}

    explicit VertexPairMap(Vertex n, T initial = T()) : n(n), values(idx({n - 2, n - 1}) + 1, initial) {}

    const_reference operator[](const VertexPair &edge) const {
        assert(edge.u < n && edge.v < n);
        return values[idx(edge)];
    }

    reference operator[](const VertexPair &edge) {
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

    [[nodiscard]] Vertex size() const { return n; }

private:
    static inline size_t idx(const VertexPair &uv) {
        return uv.v * (uv.v - 1) / 2 + uv.u;
    }
};



template <>
class VertexPairMap<bool> {
    using Block = uint_fast32_t;
    using AdjRow = boost::dynamic_bitset<Block>;
    using AdjMatrix = std::vector<AdjRow>;

    unsigned int n;
    AdjMatrix adj;

public:
    using const_reference = AdjRow::const_reference;
    class reference {
        VertexPairMap<bool>& map; VertexPair uv;
    public:
        reference(VertexPairMap<bool>& map, VertexPair uv) : map(map), uv(uv) {}
        reference& operator=(bool value) {
            map.adj[uv.u][uv.v] = true;
            map.adj[uv.v][uv.u] = true;
            return *this;
        }
        explicit operator bool() const {
            return map.adj[uv.u][uv.v];
        }
    };


    explicit VertexPairMap(const Graph &graph, bool initial = false) : VertexPairMap<bool>(graph.size()) {}

    explicit VertexPairMap(Vertex n, bool initial = false) : n(n),
                                                        adj(n, AdjRow(n)) {
        if (initial) {
            for (auto& row : adj) {
                row.set();
            }
        }
    }

    const_reference operator[](const VertexPair &uv) const {
        assert(uv.u < n && uv.v < n);
        return adj[uv.u][uv.v];
    }

    reference operator[](const VertexPair &uv) {
        assert(uv.u < n && uv.v < n);
        return reference(*this, uv);
    }

    [[nodiscard]] Vertex size() const { return n; }
};

#endif //CONCEPT_VERTEXPAIRMAP_H
