//
// Created by jonas on 27.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FIXED_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FIXED_H


#include "../Configuration.h"


/**
 * Goal: Execute fpt algorithm exactly once with editing cost determined by the configuration (k_max).
 */
class Fixed : public SearchStrategyI {
    Configuration m_config;

public:
    explicit Fixed(Configuration config) : m_config(std::move(config)) {
        if (m_config.k_max < 0) {
            throw std::runtime_error("FixedStrategy needs non-negative maximum editing cost config.k_max.");
        }
    }

    void set_initial_search_k(Cost /*initial_search_k*/) override {
        throw std::runtime_error("Fixed::set_initial_search_k not implemented.");
    }

    std::optional<Cost> get_next_search_k() override {
        return m_config.k_max;
    }

    void register_call(Cost /*k*/) override {}

    void register_bound(Cost /*k*/, Cost /*lower_bound_k*/) override {}
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FIXED_H
