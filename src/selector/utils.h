//
// Created by jonas on 04.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SELECTOR_UTILS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SELECTOR_UTILS_H


#include "FirstEditable.h"
#include "LeastWeight.h"
#include "MostMarkedPairs.h"


namespace Selector {
    using Options::Selector;

    /**
     * Factory function for Selectors. The enum selector determines the class.
     *
     * @param selector
     * @param finder
     * @param instance
     * @param marked
     * @return
     */
    static std::unique_ptr<SelectorI>
    make(Selector selector, const std::shared_ptr<FinderI> &finder, const Instance &instance,
         const VertexPairMap<bool> &marked) {
        switch (selector) {
            case Selector::LeastWeight:
                return std::make_unique<LeastWeight>(instance.costs, finder, marked);
            case Selector::FirstEditable:
                return std::make_unique<FirstEditable>(finder, marked);
            default:
                assert(false);
                return nullptr;
        }
    }
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SELECTOR_UTILS_H
