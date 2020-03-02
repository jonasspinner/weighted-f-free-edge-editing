#include <utility>

//
// Created by jonas on 02.03.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_INCREMENTBYCONSTVALUE_H
#define WEIGHTED_F_FREE_EDGE_EDITING_INCREMENTBYCONSTVALUE_H


#include "SearchStrategyI.h"


class IncrementByConstValue : public SearchStrategyI {
protected:
    Cost m_current_search_k{-1};
    Cost m_increment_value;
    Configuration m_config;

    IncrementByConstValue(Configuration config) : m_config(std::move(config)) {}

public:
    IncrementByConstValue(Cost increment_value, Configuration config) : m_increment_value(increment_value),
                                                                        m_config(std::move(config)) {}

    void set_initial_search_k(Cost initial_search_k) override {
        m_current_search_k = initial_search_k;
    }

    std::optional<Cost> get_next_search_k() override {
        Cost search_k = m_current_search_k;

        if (m_config.k_max >= 0 && search_k < m_config.k_max) {
            m_current_search_k += m_increment_value;
            return search_k;
        } else {
            return std::nullopt;
        }
    }

    void register_call(Cost /*k*/) override {}

    void register_bound(Cost /*k*/, Cost /*lower_bound_k*/) override {}
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_INCREMENTBYCONSTVALUE_H
