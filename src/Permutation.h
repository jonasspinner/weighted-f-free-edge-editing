//
// Created by jonas on 18.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_PERMUTATION_H
#define WEIGHTED_F_FREE_EDGE_EDITING_PERMUTATION_H


#include "graph/Graph.h"
#include "forbidden_subgraphs/Subgraph.h"
#include "forbidden_subgraphs/SubgraphC4P4.h"
#include "Instance.h"

class Permutation {
private:
    std::vector<Vertex> P;
    int m_seed = 0;

    explicit Permutation(std::vector<Vertex> P_, int seed) : P(std::move(P_)), m_seed(seed) {}

public:
    explicit Permutation(size_t size) : P(size) {
        for (Vertex u = 0; u < P.size(); ++u)
            P[u] = u;
    }

    Permutation(size_t size, int seed) : Permutation(size) {
        m_seed = seed;
        if (seed != 0)
            std::shuffle(P.begin(), P.end(), std::mt19937_64(static_cast<unsigned long>(seed)));
    }

    inline Vertex operator[](Vertex u) const {
        assert(u < P.size());
        return P[u];
    }

    VertexPair operator[](VertexPair uv) const {
        assert(uv.u < P.size());
        assert(uv.v < P.size());
        return {P[uv.u], P[uv.v]};
    }

    SubgraphT<Options::FSG::C4P4> operator[](const SubgraphT<Options::FSG::C4P4> &subgraph) const {
        auto result{subgraph};
        result[0] = P[result[0]];
        result[1] = P[result[1]];
        result[2] = P[result[2]];
        result[3] = P[result[3]];
        return result;
    }

    template <class T>
    std::vector<T> operator[](std::vector<T> vec) const {
        std::vector<T> result;
        result.reserve(vec.size());
        for (const T& t : vec)
            result.push_back((*this)[t]);
        return result;
    }

    template<class T>
    VertexPairMap<T> operator[](const VertexPairMap<T> &map) const {
        VertexPairMap<T> new_map(map.size());
        for (Vertex u = 0; u < map.size(); ++u) {
            for (Vertex v = u + 1; v < map.size(); ++v) {
                VertexPair uv = {u, v};
                new_map[(*this)[uv]] = map[uv];
            }
        }
#ifndef NDEBUG
        for (Vertex u = 0; u < map.size(); ++u) {
            for (Vertex v = u + 1; v < map.size(); ++v) {
                VertexPair uv = {u, v};
                assert(new_map[(*this)[uv]] == map[uv]);
            }
        }
#endif
        return new_map;
    }

    Graph operator[](const Graph &graph) const {
        Graph new_graph(graph.size());
        for (VertexPair uv : graph.edges())
            new_graph.setEdge((*this)[uv]);

#ifndef NDEBUG
        for (VertexPair uv : graph.vertexPairs())
            assert(new_graph.hasEdge((*this)[uv]) == graph.hasEdge(uv));
#endif
        return new_graph;
    }

    Instance operator[](const Instance &instance) const {
        if (instance.permutation == -m_seed)
            return Instance(instance.name, (*this)[instance.graph], (*this)[instance.costs], instance.multiplier, 0);
        else if (instance.permutation == 0)
            return Instance(instance.name, (*this)[instance.graph], (*this)[instance.costs], instance.multiplier, m_seed);
        else
            throw std::runtime_error("at most one application of permutation is allowed for instances");
    }

    [[nodiscard]] Permutation reverse() const {
        std::vector<Vertex> r_P(P.size());
        for (Vertex u = 0; u < P.size(); ++u)
            r_P[P[u]] = u;

        return Permutation(r_P, -m_seed);
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_PERMUTATION_H
