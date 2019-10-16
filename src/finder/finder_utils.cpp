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


namespace Finder {
    /**
     * Factory function for Finders. The enum forbidden determines the class.
     *
     * @param forbidden
     * @param graph
     * @return
     */
    std::unique_ptr<FinderI> make(Options::FSG forbidden, const Graph &graph) {
        using Options::FSG;
        switch (forbidden) {
            case FSG::P3:
                return std::make_unique<CenterP3>(graph);
            case FSG::C4P4:
                return std::make_unique<CenterC4P4>(graph);
            case FSG::C4_C5_2K2:
                return std::make_unique<SplitGraph>(graph);
            case FSG::C4_C5_P5_Bowtie_Necktie:
                return std::make_unique<SplitCluster>(graph);
            case FSG::C5P5:
                return std::make_unique<CenterRecC5P5>(graph);
            default:
                assert(false);
                return nullptr;
        }
    }
}
