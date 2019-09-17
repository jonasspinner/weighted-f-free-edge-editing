//
// Created by jonas on 15.07.19.
//


#include "TrivialLowerBound.h"


namespace LowerBound {
    /**
     * A trivial lower bound on the number of edits needed to solve an instance is 0.
     *
     * @return 0
     */
    Cost TrivialLowerBound::result(Cost /*k*/) {
        return 0;
    }
}
