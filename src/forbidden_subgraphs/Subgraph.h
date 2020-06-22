#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHT_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHT_H


#include "../options.h"
#include "../graph/Graph.h"
#include "../graph/VertexPairMap.h"


template<Options::FSG SetOfForbiddenSubgraphs>
class SubgraphT {
    struct SubgraphFinder{};
    class VertexIt {};
    class VertexPairIt {};
    class NonConvertingEdits{};
    using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
public:
    using Finder = SubgraphFinder;

    SubgraphT() {
        throw std::runtime_error("unimplemented");
    }

    VertexIt vertices() const;
    VertexPairIt vertex_pairs() const;
    NonConvertingEdits non_converting_edits() const;
    Cost calculate_min_cost(const VertexPairMap<Cost> &costs, const VertexPairMap<bool> &marked) const;
    [[nodiscard]] bool operator==(const Subgraph &other) const;
    [[nodiscard]] bool operator!=(const Subgraph &other) const;
    [[nodiscard]] bool operator<(const Subgraph &other) const;
    [[nodiscard]] bool contains(Vertex u) const;
    [[nodiscard]] bool contains(VertexPair uv) const;

    template<Options::FSG Set>
    friend std::ostream &operator<<(std::ostream &os, const SubgraphT<Set> &subgraph);
    template<Options::FSG Set>
    friend YAML::Emitter &operator<<(YAML::Emitter &out, const SubgraphT<Set> &subgraph);
};



#endif //WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHT_H
