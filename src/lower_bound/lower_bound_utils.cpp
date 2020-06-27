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
    std::unique_ptr<LowerBoundI>
    make(Options::LB lower_bound, const std::shared_ptr<FinderI> &finder, const Instance &instance,
         const VertexPairMap<bool> &marked, const SubgraphStats &subgraph_stats, Configuration config) {
        using Options::LB;
        switch (lower_bound) {
            case LB::Trivial:
                return std::make_unique<Trivial>();
            case LB::LocalSearch:
                return std::make_unique<LocalSearch>(instance, marked, subgraph_stats, finder);
            case LB::Greedy:
            {
                if (finder->forbidden_subgraphs() == Options::FSG::C4P4) {
                    return std::make_unique<Greedy<Options::FSG::C4P4>>(instance, marked);
                }
                throw std::runtime_error("Greedy not specialized for given forbidden subgraphs.");
            }
            case LB::SortedGreedy:
            {
                if (finder->forbidden_subgraphs() == Options::FSG::C4P4) {
                    return std::make_unique<SortedGreedy<Options::FSG::C4P4>>(instance, marked);
                }
                throw std::runtime_error("SortedGreedy not specialized for given forbidden subgraphs.");
            }
            case LB::LPRelaxation:
#ifdef GUROBI_FOUND
            {
                if (finder->forbidden_subgraphs() == Options::FSG::C4P4) {
                    return std::make_unique<LPRelaxation<Options::FSG::C4P4>>(instance, marked, std::move(config));
                }
                throw std::runtime_error("LPRelaxation not specialized for given forbidden subgraphs.");
            }
#else
                throw std::runtime_error("gurobi has to be installed to use LB::LinearProgram.");
#endif
            case LB::NPS_MWIS_Solver:
#ifdef NPS_MWIS_FOUND
            {
                if (finder->forbidden_subgraphs() == Options::FSG::C4P4) {
                    return std::make_unique<NPS_MWIS_Solver<Options::FSG::C4P4>>(instance, marked);
                }
                throw std::runtime_error("NPS_MWIS_Solver not specialized for given forbidden subgraphs.");
            }
#else
                throw std::runtime_error("nps_mwis has to be installed to use LB::NPS_MWIS_Solver.");
#endif
            case LB::LSSWZ_MWIS_Solver:
#ifdef LSSWZ_MWIS_FOUND
            {
                if (finder->forbidden_subgraphs() == Options::FSG::C4P4) {
                    return std::make_unique<LSSWZ_MWIS_Solver<Options::FSG::C4P4>>(instance, marked, std::move(config));
                }
                throw std::runtime_error("LSSWZ_MWIS_Solver not specialized for given forbidden subgraphs.");
            }
#else
                throw std::runtime_error("lsswz_mwis has to be installed to use LB::LSSWZ_MWIS_Solver.");
#endif
            case LB::GreedyWeightedPacking:
            {
                if (finder->forbidden_subgraphs() == Options::FSG::C4P4) {
                    return std::make_unique<GreedyWeightedPacking<Options::FSG::C4P4>>(instance, marked);
                }
                throw std::runtime_error("GreedyWeightedPacking not specialized for given forbidden subgraphs.");
            }
            case LB::WeightedPackingLocalSearch:
                return std::make_unique<WeightedPackingLocalSearch>(instance, marked, subgraph_stats, finder);
            default:
                throw std::runtime_error("Lower bound not found.");
        }
    }
}
