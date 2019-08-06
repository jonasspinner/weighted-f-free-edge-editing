//
// Created by jonas on 03.07.19.
//

#ifndef CONCEPT_SUBGRAPH_H
#define CONCEPT_SUBGRAPH_H


#include "Graph.h"


class Subgraph : public std::vector<Vertex> {
public:

    Subgraph(std::initializer_list<Vertex> list) : std::vector<Vertex>(list) {
#ifndef NDEBUG
        int n = 0;
        for (auto x : list)
            for (auto y : list)
                n += (x == y);
        if (n != list.size()) std::cout << *this << "\n";
        assert(n == list.size() && "vertices are not unique");
#endif
    }

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

constexpr Cost invalid_cost = std::numeric_limits<Cost>::max();

Cost cost(const Subgraph& subgraph, const VertexPairMap<bool> &marked, const VertexPairMap<Cost> &costs) {
    Cost min_cost = invalid_cost;
    subgraph.for_all_unmarked_vertex_pairs(marked, [&](VertexPair uv) {
        min_cost = std::min(min_cost, costs[uv]);
        return false;
    });
    return min_cost;
}


#endif //CONCEPT_SUBGRAPH_H
