//
// Created by jonas on 27.02.20.
//


#include "search_strategy_utils.h"

#include "Exponential.h"
#include "Fixed.h"
#include "IncrementByMinCost.h"
#include "IncrementByMultiplier.h"
#include "PrunedDelta.h"

namespace search_strategy {
    std::unique_ptr<SearchStrategyI> make(Options::FPTSearchStrategy search_strategy, Configuration config) {
        switch (search_strategy) {
            case Options::FPTSearchStrategy::Exponential:
                return std::make_unique<Exponential>();
            case Options::FPTSearchStrategy::Fixed:
                return std::make_unique<Fixed>(std::move(config));
            case Options::FPTSearchStrategy::IncrementByMinCost:
                return std::make_unique<IncrementByMinCost>();
            case Options::FPTSearchStrategy::IncrementByMultiplier:
                return std::make_unique<IncrementByMultiplier>();
            case Options::FPTSearchStrategy::PrunedDelta:
                return std::make_unique<PrunedDelta>();
            default:
                throw std::runtime_error("Options::FPTSearchStrategy::? case not handled.");
        }
    }
}
