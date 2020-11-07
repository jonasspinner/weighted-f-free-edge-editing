#include "lower_bound_utils.h"

#include "Greedy.h"
#include "GreedyWeightedPacking.h"
#include "weighted_packing/WeightedPackingLocalSearch.h"
#include "SortedGreedy.h"
#include "LocalSearch.h"
#include "Trivial.h"

#ifdef GUROBI_FOUND

#include "LPRelaxation.h"

#endif

#ifdef NPS_MWIS_FOUND

#include "NPS_MWIS_Solver.h"

#endif

#ifdef LSSWZ_MWIS_FOUND

#include "LSSWZ_MWIS_Solver.h"

#endif

namespace lower_bound {
    /**
     * A factory function for LowerBounds. The enum lower_bound determines the class.
     *
     * @param lower_bound
     * @param finder
     * @param instance
     * @param marked
     * @return
     */
    template<Options::FSG FSG>
    std::unique_ptr<LowerBoundI> make(Options::LB lower_bound, const Instance &instance,
                                      const VertexPairMap<bool> &marked, const SubgraphStats<FSG> &subgraph_stats,
                                      Configuration config) {
        using Options::LB;
        switch (lower_bound) {
            case LB::Trivial:
                return std::make_unique<Trivial>();
            case LB::LocalSearch:
                return std::make_unique<LocalSearch<FSG>>(instance, marked, subgraph_stats, config.seed);
            case LB::Greedy:
                return std::make_unique<Greedy<FSG>>(instance, marked);
            case LB::SortedGreedy:
                return std::make_unique<SortedGreedy<FSG>>(instance, marked);
            case LB::LPRelaxation:
#ifdef GUROBI_FOUND
                return std::make_unique<LPRelaxation<FSG>>(instance, marked, config.verbosity, config.timelimit);
#else
                throw std::runtime_error("gurobi has to be installed to use LB::LinearProgram.");
#endif
            case LB::NPS_MWIS_Solver:
#ifdef NPS_MWIS_FOUND
                return std::make_unique<NPS_MWIS_Solver<FSG>>(instance, marked);
#else
                throw std::runtime_error("nps_mwis has to be installed to use LB::NPS_MWIS_Solver.");
#endif
            case LB::LSSWZ_MWIS_Solver:
#ifdef LSSWZ_MWIS_FOUND
                return std::make_unique<LSSWZ_MWIS_Solver<FSG>>(instance, marked, config.verbosity);
#else
                throw std::runtime_error("lsswz_mwis has to be installed to use LB::LSSWZ_MWIS_Solver.");
#endif
            case LB::GreedyWeightedPacking:
                return std::make_unique<GreedyWeightedPacking<FSG>>(instance, marked);
            case LB::WeightedPackingLocalSearch:
                return std::make_unique<WeightedPackingLocalSearch<FSG>>(instance, marked, subgraph_stats);
            default:
                throw std::runtime_error("Lower bound not found.");
        }
    }


    template
    std::unique_ptr<LowerBoundI> make<Options::FSG::C4P4>(Options::LB lower_bound, const Instance &instance,
                                                          const VertexPairMap<bool> &marked,
                                                          const SubgraphStats<Options::FSG::C4P4> &subgraph_stats,
                                                          Configuration config);
}
