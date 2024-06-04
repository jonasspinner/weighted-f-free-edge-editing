//
// Created by jonas on 19.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FINDER_UTILS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FINDER_UTILS_H


#include "../../options.h"
#include "FinderI.h"


namespace Finder {
    std::unique_ptr<FinderI> make(Options::FSG forbidden);
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_FINDER_UTILS_H
