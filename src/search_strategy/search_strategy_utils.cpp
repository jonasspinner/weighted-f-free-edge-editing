//
// Created by jonas on 27.02.20.
//


#include "ExponentialSearchStrategy.h"
#include "FixedSearchStrategy.h"
#include "IncrementByMinCostSearchStrategy.h"
#include "IncrementByMultiplierSearchStrategy.h"
#include "PrunedDeltaSearchStrategy.h"

namespace SearchStrategy {
    std::unique_ptr<SearchStrategyI>
    make(Options::FPTSearchStrategy search_strategy);
}
