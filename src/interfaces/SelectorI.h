//
// Created by jonas on 17.07.19.
//

#ifndef CONCEPT_SELECTORI_H
#define CONCEPT_SELECTORI_H

#include "ConsumerI.h"
#include "../graph/Subgraph.h"

struct Problem {
    std::vector<VertexPair> pairs;
    bool solved = false;
};

class SelectorI : public ConsumerI {
public:
    explicit SelectorI(std::shared_ptr<FinderI> finder) : ConsumerI(std::move(finder)) {}

    virtual Problem result(StateI &state, Cost k) = 0;
};

#endif //CONCEPT_SELECTORI_H
