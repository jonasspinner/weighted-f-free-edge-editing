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
    std::unique_ptr<SelectorI> make(Options::Selector selector, const Instance &instance,
                                    const VertexPairMap<bool> &marked, const SubgraphStats<FSG> &subgraph_stats) {
        switch (selector) {
            case Options::Selector::LeastWeight:
                return std::make_unique<LeastWeight<FSG>>(instance.graph, instance.costs, marked);
            case Options::Selector::FirstFound:
                return std::make_unique<FirstFound<FSG>>(instance.graph, marked);
            case Options::Selector::MostMarkedPairs:
                return std::make_unique<MostMarkedPairs<FSG>>(instance.graph, marked, subgraph_stats);
            case Options::Selector::MostAdjacentSubgraphs:
                return std::make_unique<MostAdjacentSubgraphs<FSG>>(instance.graph, marked, subgraph_stats);
            case Options::Selector::SingleEdgeEditing:
                return std::make_unique<SingleEdgeEditing<FSG>>(instance.graph, marked, subgraph_stats);
            default:
                throw std::runtime_error("Invalid selector");
        }
    }

    template
    std::unique_ptr<SelectorI> make<Options::FSG::C4P4>(Options::Selector selector, const Instance &instance,
                                                        const VertexPairMap<bool> &marked,
                                                        const SubgraphStats<Options::FSG::C4P4> &subgraph_stats);
}
