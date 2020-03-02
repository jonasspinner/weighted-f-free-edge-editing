//
// Created by jonas on 27.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_INCREMENTBYMULTIPLIER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_INCREMENTBYMULTIPLIER_H


#include "IncrementByConstValue.h"
#include "../lower_bound/LowerBoundI.h"


class IncrementByMultiplier : public IncrementByConstValue {
public:
    IncrementByMultiplier(Configuration config) : IncrementByConstValue(std::move(config)) {
        m_increment_value = m_config.multiplier;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_INCREMENTBYMULTIPLIER_H
