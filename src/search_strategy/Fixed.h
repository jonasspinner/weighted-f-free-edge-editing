//
// Created by jonas on 27.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FIXED_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FIXED_H


#include "../Configuration.h"


class Fixed : public SearchStrategyI {
    Configuration m_config;

public:
    explicit Fixed(Configuration config) : m_config(std::move(config)) {
        if (m_config.k_max < 0) {
            throw std::runtime_error("FixedStrategy needs non-negative maximum editing cost config.k_max.");
        }
    }

    Cost search_step() override {
        return m_config.k_max;
    }

    void bound(Cost /*k*/, Cost /*lower_bound_k*/) override {}
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FIXED_H
