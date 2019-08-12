//
// Created by jonas on 29.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_GREEDYLOWERBOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_GREEDYLOWERBOUND_H


#include "../interfaces/LowerBoundI.h"
#include "../Instance.h"
#include "../graph/Subgraph.h"

namespace LowerBound {
    class GreedyLowerBound : public LowerBoundI {
        const Graph &graph;
        const VertexPairMap<Cost> &costs;
        const VertexPairMap<bool> &m_forbidden;
    public:
        class State : public StateI {
            std::unique_ptr<StateI> copy() override {
                return std::make_unique<State>(*this);
            }
        };

        GreedyLowerBound(const Instance &instance, const VertexPairMap<bool> &forbidden,
                         std::shared_ptr<FinderI> finder_ref) : LowerBoundI(std::move(finder_ref)), graph(instance.graph),
                                                            costs(instance.costs), m_forbidden(forbidden) {}

        /**
         * Calculates a lower bound on the costs required to solve the current instance.
         *
         * Greedily inserts forbidden subgraphs into the bound. Higher minimum editing costs are preferred. A subgraph
         * is not inserted if it shares an editable vertex pair with a subgraph already in the bound.
         *
         * @param state Not used
         * @param k Not used
         * @return A lower bound on the costs required to solve the current instance.
         */
        Cost result(StateI &/*state*/, Cost /*k*/) override {

            // Find all forbidden subgraphs with editable vertex pairs
            // The cost for a single forbidden subgraph is the minimum edit cost for an editable vertex pair
            std::vector<std::pair<Cost, Subgraph>> subgraphs;
            finder->find([&](Subgraph &&subgraph) {
                Cost min_cost = cost(subgraph, m_forbidden, costs);
                subgraphs.emplace_back(min_cost, std::move(subgraph));
                return false;
            });

            // Sort subgraphs with decreasing costs
            std::sort(subgraphs.begin(), subgraphs.end(),
                      [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

            if (!subgraphs.empty() && subgraphs[0].first == invalid_cost) return 0;

            Cost bound_size = 0;
            VertexPairMap<bool> is_in_bound(graph.size(), false);

            // Insert forbidden subgraphs with decreasing minimum edit cost into the bound
            // Only insert a subgraph if it does not share an editable vertex pair with a subgraph already in the bound
            for (const auto&[cost, subgraph] : subgraphs) {
                bool touches_bound = false;
                for (VertexPair uv : subgraph.vertexPairs()) {
                    if (!m_forbidden[uv] && is_in_bound[uv]) touches_bound = true;
                }

                if (!touches_bound) {
                    bound_size += cost;
                    for (VertexPair uv : subgraph.vertexPairs()) {
                        if (!m_forbidden[uv]) is_in_bound[uv] = true;
                    }
                }
            }

            return bound_size;
        }

        std::unique_ptr<StateI> initialize(Cost /*k*/) override { return std::make_unique<State>(); }

        void before_mark_and_edit(StateI &/*state*/, VertexPair /*uv*/) override {}

        void after_mark_and_edit(StateI &/*state*/, VertexPair /*uv*/) override {}

        void before_mark(StateI &/*state*/, VertexPair /*uv*/) override {}

        void after_mark(StateI &/*state*/, VertexPair /*uv*/) override {}

        void before_edit(StateI &/*state*/, VertexPair /*uv*/) override {}

        void after_edit(StateI &/*state*/, VertexPair /*uv*/) override {}

        void after_unmark(StateI &/*state*/, VertexPair /*uv*/) override {}
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_GREEDYLOWERBOUND_H
