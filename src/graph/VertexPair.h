//
// Created by jonas on 12.03.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIR_H
#define WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIR_H


#include "../definitions.h"


class VertexPair {
public:
    Vertex u;
    Vertex v;

    /**
     * An unordered pair of vertices. Represents an undirected edge.
     *
     * The vertex pairs {u, v} and {v, u} are the same object. It is guaranteed that pair.u < pair.v.
     *
     * @param x A vertex
     * @param y A vertex
     */
    VertexPair(Vertex x, Vertex y) : u(x < y ? x : y), v(x < y ? y : x) { assert(x != y); }

    struct Ordered {};
    VertexPair(Vertex x, Vertex y, Ordered) : u(x), v(y) { assert(x < y); } // x < y as precondition

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


#endif //WEIGHTED_F_FREE_EDGE_EDITING_VERTEXPAIR_H
