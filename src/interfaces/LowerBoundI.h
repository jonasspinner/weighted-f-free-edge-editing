//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_LOWERBOUNDI_H
#define CONCEPT_LOWERBOUNDI_H


#include "ConsumerI.h"


class LowerBoundI : public ConsumerI {
public:
    explicit LowerBoundI(std::shared_ptr<FinderI> finder) : ConsumerI(std::move(finder)) {}

    virtual Cost result(StateI &state, Cost k) = 0;
};

#endif //CONCEPT_LOWERBOUNDI_H
