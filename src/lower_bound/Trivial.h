//
// Created by jonas on 15.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_TRIVIAL_H
#define WEIGHTED_F_FREE_EDGE_EDITING_TRIVIAL_H


#include "LowerBoundI.h"


namespace lower_bound {
    class Trivial : public LowerBoundI {
    public:
        explicit Trivial(std::shared_ptr<FinderI> finder_ptr) : LowerBoundI(std::move(finder_ptr)) {};

        Cost calculate_lower_bound(Cost /*k*/) override;
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_TRIVIAL_H
