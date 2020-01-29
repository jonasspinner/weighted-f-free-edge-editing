//
// Created by jonas on 06.09.19.
//


#include "lower_bound_utils.h"

#include "GreedyLowerBound.h"
#include "SortedGreedyLowerBound.h"
#include "LocalSearchLowerBound.h"
#include "TrivialLowerBound.h"

#ifdef GUROBI_FOUND

#include "LPRelaxationLowerBound.h"

#endif

#include "ILSMWISSolverLowerBound.h"

namespace LowerBound {
    /**
     * A factory function for LowerBounds. The enum lower_bound determines the class.
     *
     * @param lower_bound
     * @param finder
     * @param instance
     * @param marked
     * @return
     */
    std::unique_ptr<LowerBoundI>
    make(Options::LB lower_bound, const std::shared_ptr<FinderI> &finder, const Instance &instance,
         const VertexPairMap<bool> &marked, const SubgraphStats &subgraph_stats, Configuration config) {
        using Options::LB;
        switch (lower_bound) {
            case LB::Trivial:
                return std::make_unique<TrivialLowerBound>(finder);
            case LB::LocalSearch:
                return std::make_unique<LocalSearchLowerBound>(instance, marked, subgraph_stats, finder);
            case LB::Greedy:
                return std::make_unique<GreedyLowerBound>(instance, marked, finder);
            case LB::SortedGreedy:
                return std::make_unique<SortedGreedyLowerBound>(instance, marked, finder);
            case LB::LPRelaxation:
#ifdef GUROBI_FOUND
            return std::make_unique<LPRelaxationLowerBound>(instance, marked, std::move(config), finder);
#else
            std::cerr << "gurobi has to be installed to use LB::LinearProgram.";
                abort();
#endif
            case LB::ILSMWISSolver:
                return std::make_unique<ILSMWISSolverLowerBound>(instance, marked, finder);
            default:
                assert(false);
                return nullptr;
        }
    }
}
