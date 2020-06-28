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
    std::unique_ptr<SelectorI>
    make(Options::Selector selector, Options::FSG fsg, const Instance &instance,
         const VertexPairMap<bool> &marked, const SubgraphStats &subgraph_stats) {
        switch (selector) {
            case Options::Selector::LeastWeight: {
                switch (fsg) {
                    case Options::FSG::C4P4:
                        return std::make_unique<LeastWeight<Options::FSG::C4P4>>(instance.graph, instance.costs, marked);
                    default:
                        throw std::runtime_error("LeastWeight not specialized for given forbidden subgraphs.");
                }
            }
            case Options::Selector::FirstFound: {
                switch (fsg) {
                    case Options::FSG::C4P4:
                        return std::make_unique<FirstFound<Options::FSG::C4P4>>(instance.graph, marked);
                    default:
                        throw std::runtime_error("FirstFound not specialized for given forbidden subgraphs.");
                }
            }
            case Options::Selector::MostMarkedPairs: {
                switch (fsg) {
                    case Options::FSG::C4P4:
                    return std::make_unique<MostMarkedPairs<Options::FSG::C4P4>>(
                            instance.graph, marked, subgraph_stats);
                    default:
                        throw std::runtime_error("MostMarkedPairs not specialized for given forbidden subgraphs.");
                }
            }
            case Options::Selector::MostAdjacentSubgraphs: {
                switch (fsg) {
                    case Options::FSG::C4P4:
                        return std::make_unique<MostAdjacentSubgraphs<Options::FSG::C4P4>>(
                                instance.graph, marked, subgraph_stats);
                    default:
                        throw std::runtime_error("MostAdjacentSubgraphs not specialized for given forbidden subgraphs.");
                }
            }
            case Options::Selector::SingleEdgeEditing: {
                switch (fsg) {
                    case Options::FSG::C4P4:
                    return std::make_unique<SingleEdgeEditing<Options::FSG::C4P4>>(
                            instance.graph, marked, subgraph_stats);
                    default:
                        throw std::runtime_error("SingleEdgeEditing not specialized for given forbidden subgraphs.");
                }
            }
            default:
                throw std::runtime_error("Invalid selector");
        }
    }
}
