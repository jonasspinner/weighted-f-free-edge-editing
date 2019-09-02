//
// Created by jonas on 15.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_TRIVIALLOWERBOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_TRIVIALLOWERBOUND_H


#include "../interfaces/LowerBoundI.h"


namespace LowerBound {
    class TrivialLowerBound : public LowerBoundI {
    public:
        explicit TrivialLowerBound(std::shared_ptr<FinderI> finder_ptr) : LowerBoundI(std::move(finder_ptr)) {};

        Cost result(Cost /*k*/) override {
            return 0;
        }
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_TRIVIALLOWERBOUND_H
