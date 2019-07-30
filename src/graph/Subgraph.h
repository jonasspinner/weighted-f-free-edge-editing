//
// Created by jonas on 03.07.19.
//

#ifndef CONCEPT_SUBGRAPH_H
#define CONCEPT_SUBGRAPH_H


#include <array>

#include "Graph.h"


class Subgraph : public std::vector<Vertex> {
public:

    Subgraph(std::initializer_list<Vertex> list) : std::vector<Vertex>(list) {}

    template <typename VertexPairCallback>
    bool for_all_vertex_pairs(VertexPairCallback callback) const {
        for (size_t i = 0; i < size(); ++i) {
            for (size_t j = i+1; j < size(); ++j) {
                if(callback(VertexPair((*this)[i], (*this)[j]))) return true;
            }
        }
        return false;
    }

    template <typename VertexPairCallback>
    bool for_all_unmarked_vertex_pairs(const VertexPairMap<bool> &marked, VertexPairCallback callback) const {
        return for_all_vertex_pairs([&](VertexPair uv) {
            return !marked[uv] && callback(uv);
        });
    }

    friend std::ostream&operator<<(std::ostream& os, const Subgraph& subgraph) {
        os << "{";
        for (Vertex u : subgraph) os << " " << u;
        return os << " }";
    }
};


#endif //CONCEPT_SUBGRAPH_H
