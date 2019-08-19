//
// Created by jonas on 18.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_PERMUTATION_H
#define WEIGHTED_F_FREE_EDGE_EDITING_PERMUTATION_H


class Permutation {
private:
    std::vector<Vertex> P;

    explicit Permutation(std::vector<Vertex> P_) : P(std::move(P_)) {}

public:
    explicit Permutation(size_t size) : P(size) {
        for (Vertex u = 0; u < P.size(); ++u)
            P[u] = u;
    }

    Permutation(size_t size, int seed) : Permutation(size) {
        std::shuffle(P.begin(), P.end(), std::mt19937_64(seed));
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

        return new_map;
    }

    Graph operator[](const Graph &graph) const {
        Graph new_graph(graph.size());
        for (VertexPair uv : graph.edges())
            new_graph.set_edge((*this)[uv]);

        return new_graph;
    }

    Instance operator[](const Instance &instance) const {
        return {(*this)[instance.graph], (*this)[instance.costs]};
    }

    [[nodiscard]] Permutation reverse() const {
        std::vector<Vertex> r_P(P.size());
        for (Vertex u = 0; u < P.size(); ++u)
            r_P[P[u]] = u;

        return Permutation(r_P);
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_PERMUTATION_H
