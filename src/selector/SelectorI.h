//
// Created by jonas on 17.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SELECTORI_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SELECTORI_H

#include "../consumer/ConsumerI.h"
#include "../graph/Subgraph.h"


struct Problem {
    std::vector<VertexPair> pairs;
    bool solved = false;
};

class SelectorI : public ConsumerI {
protected:
    std::shared_ptr<FinderI> finder;
public:
    explicit SelectorI(std::shared_ptr<FinderI> finder_ptr) : finder(std::move(finder_ptr)) {}

    virtual Problem select_problem(Cost k) = 0;
};

#endif //WEIGHTED_F_FREE_EDGE_EDITING_SELECTORI_H
