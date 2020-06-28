#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SELECTOR_UTILS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SELECTOR_UTILS_H


#include <memory>

#include "SelectorI.h"
#include "../options.h"
#include "../Instance.h"
#include "../consumer/SubgraphStats.h"


namespace selector {
    std::unique_ptr<SelectorI>
    make(Options::Selector selector, Options::FSG fsg, const Instance &instance,
                   const VertexPairMap<bool> &marked, const SubgraphStats &subgraph_stats);

}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_SELECTOR_UTILS_H
