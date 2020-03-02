//
// Created by jonas on 27.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_INCREMENTBYMINCOST_H
#define WEIGHTED_F_FREE_EDGE_EDITING_INCREMENTBYMINCOST_H


#include "IncrementByConstValue.h"


class IncrementByMinCost : public IncrementByConstValue {
public:
    IncrementByMinCost(const Instance &instance, Configuration config) : IncrementByConstValue(std::move(config)) {
        Cost min_cost = std::numeric_limits<Cost>::max();
        for (auto uv : instance.graph.vertexPairs()) {
            auto cost = instance.costs[uv];
            if (cost != 0 && cost < min_cost)
                min_cost = cost;
        }
        m_increment_value = std::max(min_cost, 1);
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_INCREMENTBYMINCOST_H
