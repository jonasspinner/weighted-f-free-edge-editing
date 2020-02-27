//
// Created by jonas on 27.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_EXPONENTIAL_H
#define WEIGHTED_F_FREE_EDGE_EDITING_EXPONENTIAL_H

#include <deque>
#include "SearchStrategyI.h"

class Exponential : public SearchStrategyI {
    std::deque<size_t> m_calls_history;

public:
    Cost search_step() override { return -1; }
    void bound(Cost /*k*/, Cost /*lower_bound_k*/) override {}
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_EXPONENTIAL_H
