//
// Created by jonas on 04.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LB_UTILS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LB_UTILS_H


#include "GreedyLowerBound.h"
#include "LocalSearchLowerBound.h"
#include "TrivialLowerBound.h"
#ifdef GUROBI_FOUND
#include "LinearProgramLowerBound.h"
#endif



namespace LowerBound {
    using Options::LB;
    
    /**
     * A factory function for LowerBounds. The enum lower_bound determines the class.
     * 
     * @param lower_bound 
     * @param finder 
     * @param instance 
     * @param marked 
     * @return 
     */
    static std::unique_ptr<LowerBoundI> make(LB lower_bound, const std::shared_ptr<FinderI> &finder, const Instance &instance, const VertexPairMap<bool> &marked) {
        switch (lower_bound) {
            case LB::No:
                return std::make_unique<TrivialLowerBound>(finder);
            case LB::LocalSearch:
                return std::make_unique<LocalSearchLowerBound>(instance, marked, finder);
            case LB::Greedy:
                return std::make_unique<GreedyLowerBound>(instance, marked, finder);
            case LB::LinearProgram:
#ifdef GUROBI_FOUND
                return std::make_unique<LinearProgramLowerBound>(instance, marked, finder);
#else
                std::cerr << "gurobi has to be installed to use LB::LinearProgram.";
                abort();
#endif
            default:
                assert(false);
                return nullptr;
        }
    }
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_LB_UTILS_H
