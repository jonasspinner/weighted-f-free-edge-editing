//
// Created by jonas on 19.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SOLVER_H


#include "../Instance.h"
#include "../Solution.h"


class Solver {
public:
    class Result {
    public:
        enum class Status {
            Solved, Unsolved, Timeout
        } status;
        std::vector<Solution> solutions;
        Result(Status status_, std::vector<Solution> solutions_) : status(status_), solutions(std::move(solutions_)) {}
        static Result Solved(std::vector<Solution> solutions_) { return Result(Status::Solved, std::move(solutions_)); }
        static Result Unsolved() { return Result(Status::Unsolved, {}); }
        static Result Timeout() { return Result(Status::Timeout, {}); }
    };

    virtual ~Solver() = default;

    virtual Result solve(Instance instance) = 0;
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SOLVER_H
