//
// Created by jonas on 06.09.19.
//


#include "finder_utils.h"

#include "Center.h"
#include "CenterC4P4.h"
#include "CenterP3.h"
#include "NaiveC4P4.h"
#include "NaiveP3.h"
#include "SplitCluster.h"
#include "SplitGraph.h"
#include "OuterP3.h"


namespace Finder {
    /**
     * Factory function for Finders. The enum forbidden determines the class.
     *
     * @param forbidden
     * @return
     */
    std::unique_ptr<FinderI> make(Options::FSG forbidden) {
        using Options::FSG;
        switch (forbidden) {
            case FSG::P3:
                return std::make_unique<OuterP3>();
            case FSG::P4:
                return std::make_unique<CenterRecP4>();
            case FSG::P5:
                return std::make_unique<CenterRecP5>();
            case FSG::P6:
                return std::make_unique<CenterRecP6>();
            case FSG::C4P4:
                return std::make_unique<CenterC4P4>();
            case FSG::C5P5:
                return std::make_unique<CenterRecC5P5>();
            case FSG::C6P6:
                return std::make_unique<CenterRecC6P6>();
            case FSG::C4_C5_2K2:
                return std::make_unique<SplitGraph>();
            case FSG::C4_C5_P5_Bowtie_Necktie:
                return std::make_unique<SplitCluster>();
            default:
                throw std::runtime_error("Finder not supported");
                return nullptr;
        }
    }
}
