#include "selector_utils.h"

#include "FirstFound.h"
#include "LeastWeight.h"
#include "MostMarkedPairs.h"
#include "MostAdjacentSubgraphs.h"
#include "SingleEdgeEditing.h"


namespace selector {

    /**
     * Factory function for Selectors. The enum selector determines the class.
     *
     * @param selector
     * @param finder
     * @param instance
     * @param marked
     * @return
     */
    template<Options::FSG FSG>
    std::unique_ptr<SelectorI> make(Options::Selector selector, const EditState* edit_state, const SubgraphStats<FSG> *subgraph_stats) {
        switch (selector) {
            case Options::Selector::LeastWeight:
                return std::make_unique<LeastWeight<FSG>>(edit_state);
            case Options::Selector::FirstFound:
                return std::make_unique<FirstFound<FSG>>(edit_state);
            case Options::Selector::MostMarkedPairs:
                return std::make_unique<MostMarkedPairs<FSG>>(edit_state, subgraph_stats);
            case Options::Selector::MostAdjacentSubgraphs:
                return std::make_unique<MostAdjacentSubgraphs<FSG>>(edit_state, subgraph_stats);
            case Options::Selector::SingleEdgeEditing:
                return std::make_unique<SingleEdgeEditing<FSG>>(edit_state, subgraph_stats);
            default:
                throw std::runtime_error("Invalid selector");
        }
    }

    template
    std::unique_ptr<SelectorI> make<Options::FSG::C4P4>(Options::Selector selector, const EditState *edit_state,
                                                        const SubgraphStats<Options::FSG::C4P4> *subgraph_stats);
}
