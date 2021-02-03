#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHT_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHT_H


#include "../options.h"
#include "../graph/Graph.h"
#include "../graph/VertexPairMap.h"


template<Options::FSG SetOfForbiddenSubgraphs>
class SubgraphT;


namespace subgraph_iterators {
    enum class IterationControl : bool {
        Continue = 0,
        Break = 1,
    };

    constexpr IterationControl break_if(bool condition) noexcept {
        if (condition)
            return IterationControl::Break;
        return IterationControl::Continue;
    }

    enum class IterationExit : bool {
        Normal = 0,
        Break = 1,
    };
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_SUBGRAPHT_H
