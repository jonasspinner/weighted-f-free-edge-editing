//
// Created by jonas on 15.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LOWERBOUNDI_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LOWERBOUNDI_H


#include "../consumer/ConsumerI.h"


class LowerBoundI : public ConsumerI {
public:
    explicit LowerBoundI(std::shared_ptr<FinderI> finder_ptr) : ConsumerI(std::move(finder_ptr)) {}

    /**
     * Calculate a lower bound for the current state. The search can be stopped when the lower bound is larger than k.
     *
     * @param k The remaining editing cost
     * @return
     */
    virtual Cost calculate_lower_bound(Cost k) = 0;

    /**
     * Returns a lower bound for the current state. The call should not improve the lower bound and must have O(1)
     * complexity.
     *
     * @return
     */
    virtual Cost get_lower_bound() { return 0; }
};

#endif //WEIGHTED_F_FREE_EDGE_EDITING_LOWERBOUNDI_H
