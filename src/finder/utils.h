//
// Created by jonas on 19.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_UTIL_H
#define WEIGHTED_F_FREE_EDGE_EDITING_UTIL_H


#include "NaiveP3.h"
#include "CenterC4P4.h"


namespace Finder {
    static std::shared_ptr<FinderI> make(Options::FSG forbidden, const Graph &graph) {
        switch (forbidden) {
            case Options::FSG::P3:
                return std::make_shared<Finder::CenterP3>(graph);
            case Options::FSG::P4C4:
                return std::make_shared<Finder::CenterC4P4>(graph);
            case Options::FSG::C4_C5_2K2:
                return std::make_shared<Finder::SplitGraph>(graph);
            case Options::FSG::C4_C5_P5_Bowtie_Necktie:
                return std::make_shared<Finder::SplitCluster>(graph);
            default:
                assert(false);
                return nullptr;
        }
    }
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_UTIL_H
