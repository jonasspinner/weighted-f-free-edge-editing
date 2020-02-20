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

#ifdef NPS_MWIS_FOUND
#include "NPS_MWIS_SolverLowerBound.h"
#endif

#ifdef LSSWZ_MWIS_FOUND
#include "LSSWZ_MWIS_SolverLowerBound.h"
#endif

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
                throw std::runtime_error("gurobi has to be installed to use LB::LinearProgram.");
#endif
            case LB::NPS_MWIS_Solver:
#ifdef NPS_MWIS_FOUND
                return std::make_unique<NPS_MWIS_SolverLowerBound>(instance, marked, finder);
#else
                throw std::runtime_error("nps_mwis has to be installed to use LB::NPS_MWIS_Solver.");
#endif
            case LB::LSSWZ_MWIS_Solver:
#ifdef LSSWZ_MWIS_FOUND
                return std::make_unique<LSSWZ_MWIS_SolverLowerBound>(instance, marked, std::move(config), finder);
#else
                throw std::runtime_error("lsswz_mwis has to be installed to use LB::LSSWZ_MWIS_Solver.");
#endif
            default:
                assert(false);
                return nullptr;
        }
    }
}
