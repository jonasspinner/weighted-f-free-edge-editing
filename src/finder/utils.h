//
// Created by jonas on 19.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FINDER_UTILS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FINDER_UTILS_H


#include "Center.h"
#include "CenterC4P4.h"
#include "CenterP3.h"
#include "NaiveC4P4.h"
#include "NaiveP3.h"
#include "SplitCluster.h"
#include "SplitGraph.h"


namespace Finder {
    using Options::FSG;

    /**
     * Factory function for Finders. The enum forbidden determines the class.
     * 
     * @param forbidden 
     * @param graph 
     * @return 
     */
    static std::unique_ptr<FinderI> make(FSG forbidden, const Graph &graph) {
        switch (forbidden) {
            case FSG::P3:
                return std::make_unique<CenterP3>(graph);
            case FSG::P4C4:
                return std::make_unique<CenterC4P4>(graph);
            case FSG::C4_C5_2K2:
                return std::make_unique<SplitGraph>(graph);
            case FSG::C4_C5_P5_Bowtie_Necktie:
                return std::make_unique<SplitCluster>(graph);
            default:
                assert(false);
                return nullptr;
        }
    }
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_FINDER_UTILS_H
