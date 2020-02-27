//
// Created by jonas on 15.07.19.
//


#include "TrivialLowerBound.h"


namespace lower_bound {
    /**
     * A trivial lower bound on the number of edits needed to solve an instance is 0.
     *
     * @return 0
     */
    Cost TrivialLowerBound::calculate_lower_bound(Cost /*k*/) {
        return 0;
    }
}
