#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LOWERBOUNDI_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LOWERBOUNDI_H


#include "../consumer/ConsumerI.h"


class LowerBoundI : public ConsumerI {
protected:
    std::shared_ptr<FinderI> finder;
public:
    explicit LowerBoundI(std::shared_ptr<FinderI> finder_ptr) : finder(std::move(finder_ptr)) {}

    /**
     * Calculate a lower bound for the current state. The search can be stopped when the lower bound is larger than k.
     *
     * @param k The remaining editing cost
     * @return
     */
    virtual Cost calculate_lower_bound(Cost k) = 0;

    virtual Cost calculate_lower_bound_no_edit_branch() {
        return 0;
    }
};

#endif //WEIGHTED_F_FREE_EDGE_EDITING_LOWERBOUNDI_H
