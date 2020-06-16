#ifndef WEIGHTED_F_FREE_EDGE_EDITING_TRIVIAL_H
#define WEIGHTED_F_FREE_EDGE_EDITING_TRIVIAL_H


#include "LowerBoundI.h"


namespace lower_bound {
    class Trivial : public LowerBoundI {
    public:
        /**
         * A trivial lower bound on the number of edits needed to solve an instance is 0.
         *
         * @return 0
         */
        Cost calculate_lower_bound(Cost /*k*/) override {
            return 0;
        }
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_TRIVIAL_H
