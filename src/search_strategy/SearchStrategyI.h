//
// Created by jonas on 27.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SEARCHSTRATEGYI_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SEARCHSTRATEGYI_H


#include <optional>
#include "../definitions.h"


/**
 * Concept
 * - A search strategy handles the search steps of the fpt algorithm.
 * - Problem: Strategies use information from previous steps.
 * - Problem: A search strategies must be able to disallow the initial lower bound.
 *      # TODO: currently not covered by design
 * - Maybe give strategy the instance and a lower bound algorithm to decide the first search step.
 */
class SearchStrategyI {
public:
    virtual ~SearchStrategyI() = default;

    virtual void set_initial_search_k(Cost search_k) = 0;
    virtual std::optional<Cost> get_next_search_k() = 0;

    virtual void register_call(Cost k) = 0;
    virtual void register_bound(Cost k, Cost lower_bound_k) = 0;
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SEARCHSTRATEGYI_H
