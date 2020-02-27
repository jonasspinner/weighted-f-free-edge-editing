//
// Created by jonas on 27.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_EXPONENTIALSEARCHSTRATEGY_H
#define WEIGHTED_F_FREE_EDGE_EDITING_EXPONENTIALSEARCHSTRATEGY_H

#include <deque>
#include "../interfaces/SearchStrategyI.h"

class ExponentialSearchStrategy : public SearchStrategyI {
    std::deque<size_t> m_calls_history;
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_EXPONENTIALSEARCHSTRATEGY_H
