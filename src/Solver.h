//
// Created by jonas on 19.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SOLVER_H


#include "Instance.h"
#include "Solution.h"


class Solver {
public:
    virtual ~Solver() = default;

    virtual std::vector<Solution> solve(Instance instance) = 0;
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SOLVER_H
