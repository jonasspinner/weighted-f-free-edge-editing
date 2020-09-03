#ifndef WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIR_H
#define WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIR_H


#include <cassert>

#include <yaml-cpp/yaml.h>

#include "../definitions.h"


class VertexPair {
public:
    Vertex u;
    Vertex v;

    constexpr VertexPair() noexcept: u(0), v(1) {}

    /**
     * An unordered pair of vertices. Represents an undirected edge.
     *
     * The vertex pairs {u, v} and {v, u} are the same object. It is guaranteed that pair.u < pair.v.
     *
     * @param x A vertex
     * @param y A vertex
     */
    constexpr VertexPair(Vertex x, Vertex y) noexcept: u(x < y ? x : y), v(x < y ? y : x) { assert(x != y); }

    struct Ordered {};
    constexpr VertexPair(Vertex x, Vertex y, Ordered) noexcept: u(x), v(y) { assert(x < y); } // x < y as precondition

    friend std::ostream &operator<<(std::ostream &os, VertexPair uv) {
        return os << "{" << uv.u << ", " << uv.v << "}";
    }

    friend YAML::Emitter &operator<<(YAML::Emitter &out, VertexPair uv) {
        out << YAML::Flow;
        out << YAML::BeginSeq << uv.u << uv.v << YAML::EndSeq;
        return out;
    }

    [[nodiscard]] constexpr friend bool operator==(const VertexPair &uv, const VertexPair &xy) {
        return std::tie(uv.u, uv.v) == std::tie(xy.u, xy.v);
    }

    [[nodiscard]] constexpr friend bool operator!=(const VertexPair &uv, const VertexPair &xy) {
        return std::tie(uv.u, uv.v) != std::tie(xy.u, xy.v);
    }

    [[nodiscard]] constexpr friend bool operator<(const VertexPair &uv, const VertexPair &xy) {
        return std::tie(uv.u, uv.v) < std::tie(xy.u, xy.v);
    }

    [[nodiscard]] constexpr friend bool operator<=(const VertexPair &uv, const VertexPair &xy) {
        return std::tie(uv.u, uv.v) <= std::tie(xy.u, xy.v);
    }

    [[nodiscard]] constexpr friend bool operator>(const VertexPair &uv, const VertexPair &xy) {
        return std::tie(uv.u, uv.v) > std::tie(xy.u, xy.v);
    }

    [[nodiscard]] constexpr friend bool operator>=(const VertexPair &uv, const VertexPair &xy) {
        return std::tie(uv.u, uv.v) >= std::tie(xy.u, xy.v);
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIR_H
