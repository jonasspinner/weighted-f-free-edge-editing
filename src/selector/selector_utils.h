#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SELECTOR_UTILS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SELECTOR_UTILS_H


#include <memory>

#include "SelectorI.h"
#include "../options.h"
#include "../Instance.h"
#include "../consumer/SubgraphStats.h"


namespace selector {
    template<Options::FSG FSG>
    std::unique_ptr<SelectorI> make(Options::Selector selector, const EditState *edit_state,
                                    const SubgraphStats<FSG> *subgraph_stats);

}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_SELECTOR_UTILS_H
