//
// Created by jonas on 27.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SEARCH_STRATEGY_UTILS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SEARCH_STRATEGY_UTILS_H


#include <memory>
#include "SearchStrategyI.h"
#include "../options.h"

namespace search_strategy {
    std::unique_ptr<SearchStrategyI> make(Options::FPTSearchStrategy search_strategy);
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_SEARCH_STRATEGY_UTILS_H
