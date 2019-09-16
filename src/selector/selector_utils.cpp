//
// Created by jonas on 06.09.19.
//


#include "selector_utils.h"

#include "FirstEditable.h"
#include "LeastWeight.h"
#include "MostMarkedPairs.h"
#include "../Configuration.h"
#include "../Instance.h"


namespace Selector {

    /**
     * Factory function for Selectors. The enum selector determines the class.
     *
     * @param selector
     * @param finder
     * @param instance
     * @param marked
     * @return
     */
    std::unique_ptr<SelectorI>
    make(Options::Selector selector, const std::shared_ptr<FinderI> &finder, const Instance &instance,
         const VertexPairMap<bool> &marked) {
        switch (selector) {
            case Options::Selector::LeastWeight:
                return std::make_unique<LeastWeight>(instance.costs, finder, marked);
            case Options::Selector::FirstEditable:
                return std::make_unique<FirstEditable>(finder, marked);
            case Options::Selector::MostMarked:
                return std::make_unique<MostMarkedPairs>(finder, marked);
            default:
                assert(false);
                return nullptr;
        }
    }
}
