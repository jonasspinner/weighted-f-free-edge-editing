//
// Created by jonas on 19.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_UTIL_H
#define WEIGHTED_F_FREE_EDGE_EDITING_UTIL_H


#include "NaiveP3.h"
#include "CenterC4P4.h"


namespace Finder {
    static std::unique_ptr<FinderI> make(Options::FSG fsg, const Graph &graph) {
        switch (fsg) {
            case Options::FSG::P3:
                return std::make_unique<Finder::NaiveP3>(graph);
            case Options::FSG::P4C4:
                return std::make_unique<Finder::CenterC4P4>(graph);
            default:
                assert(false);
                return nullptr;
        }
    }
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_UTIL_H
