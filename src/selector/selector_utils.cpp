//
// Created by jonas on 06.09.19.
//


#include "selector_utils.h"

#include "FirstFound.h"
#include "LeastWeight.h"
#include "MostMarkedPairs.h"
#include "MostAdjacentSubgraphs.h"
#include "Most.h"

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
         const VertexPairMap<bool> &marked, const SubgraphStats &subgraph_stats) {
        switch (selector) {
            case Options::Selector::LeastWeight:
                return std::make_unique<LeastWeight>(instance.costs, finder, marked);
            case Options::Selector::FirstFound:
                return std::make_unique<FirstFound>(finder, marked);
            case Options::Selector::MostMarkedPairs:
                return std::make_unique<MostMarkedPairs>(finder, marked, subgraph_stats);
            case Options::Selector::MostAdjacentSubgraphs:
                return std::make_unique<Most>(finder, marked, subgraph_stats);
            default:
                assert(false);
                return nullptr;
        }
    }
}
