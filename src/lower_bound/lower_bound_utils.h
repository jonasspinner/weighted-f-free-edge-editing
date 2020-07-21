#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LB_UTILS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LB_UTILS_H


#include "../options.h"
#include "LowerBoundI.h"
#include "../Instance.h"
#include "../consumer/SubgraphStats.h"
#include "../finder/FinderI.h"


namespace lower_bound {
    template<Options::FSG FSG>
    std::unique_ptr<LowerBoundI> make(Options::LB lower_bound, const Instance &instance,
         const VertexPairMap<bool> &marked, const SubgraphStats<FSG> &subgraph_stats, Configuration config);
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_LB_UTILS_H
