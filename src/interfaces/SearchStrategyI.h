//
// Created by jonas on 27.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SEARCHSTRATEGYI_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SEARCHSTRATEGYI_H


#include "../definitions.h"

class SearchStrategyI {
public:
    virtual ~SearchStrategyI() = default;

    virtual Cost search_step() = 0;
    virtual void bound(Cost k, Cost lower_bound_k) = 0;
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SEARCHSTRATEGYI_H
