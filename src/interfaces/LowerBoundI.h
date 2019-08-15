//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_LOWERBOUNDI_H
#define CONCEPT_LOWERBOUNDI_H


#include "ConsumerI.h"


class LowerBoundI : public ConsumerI {
public:
    explicit LowerBoundI(std::shared_ptr<FinderI> finder_ptr) : ConsumerI(std::move(finder_ptr)) {}

    virtual Cost result(Cost k) = 0;
};

#endif //CONCEPT_LOWERBOUNDI_H
