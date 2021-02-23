#include "SortedGreedy.h"

#include <vector>


namespace lower_bound {
    /**
     * Calculates a lower bound on the costs required to solve the current instance.
     *
     * Greedily inserts forbidden subgraphs into the bound. Higher minimum editing costs are preferred. A subgraph
     * is not inserted if it shares an editable vertex pair with a subgraph already in the bound.
     *
     * Complexity:
     *      m = #forbidden subgraphs
     *      O(m * log m) time
     *      O(m) additional space
     *
     * @param k Not used
     * @return A lower bound on the costs required to solve the current instance.
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    Cost SortedGreedy<SetOfForbiddenSubgraphs>::calculate_lower_bound(Cost k) {
        Cost max_min_cost = std::numeric_limits<Cost>::min();

        // Find all forbidden subgraphs with editable vertex pairs.
        // The cost for a single forbidden subgraph is the minimum edit cost for an editable vertex pair.
        std::vector<std::pair<Cost, Subgraph>> subgraphs;

        m_finder.find(m_edit_state->graph(), [&](Subgraph subgraph) {
            Cost min_cost = subgraph.calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map());
            subgraphs.emplace_back(min_cost, std::move(subgraph));
            max_min_cost = std::max(max_min_cost, min_cost);
            return subgraph_iterators::break_if(max_min_cost > k);
        });

        // If one subgraph has an already large enough cost or if a subgraph is fully marked.
        if (max_min_cost > k)
            return max_min_cost;

        // Sort subgraphs with decreasing costs.
        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

        Cost bound_size = 0;
        m_used_in_bound.clear();

        // Insert forbidden subgraphs with decreasing minimum edit cost into the bound.
        // Only insert a subgraph if it does not share an editable vertex pair with a subgraph already in the bound.
        for (const auto&[cost, subgraph] : subgraphs) {
            if (bound_size > k)
                break;

            // Check if the subgraph is adjacent to one already used in the bound.
            bool touches_bound = false;
            for (auto uv : subgraph.non_converting_edits()) {
                if (!m_edit_state->is_marked(uv) && m_used_in_bound[uv]) {
                    touches_bound = true;
                    break;
                }
            }

            if (!touches_bound) {
                bound_size += cost;
                for (auto uv : subgraph.non_converting_edits()) {
                    if (!m_edit_state->is_marked(uv))
                        m_used_in_bound[uv] = true;
                }
            }
        }

        return bound_size;
    }

    template<Options::FSG SetOfForbiddenSubgraphs>
    [[nodiscard]] std::tuple<Cost, VertexPairMap<bool>, std::vector<typename SortedGreedy<SetOfForbiddenSubgraphs>::Subgraph>, std::vector<VertexPair>>
    SortedGreedy<SetOfForbiddenSubgraphs>::calculate_lower_bound_and_packing() {

        // Find all forbidden subgraphs with editable vertex pairs.
        // The cost for a single forbidden subgraph is the minimum edit cost for an editable vertex pair.
        std::vector<std::pair<Cost, Subgraph>> subgraphs;

        m_finder.find(m_edit_state->graph(), [&](Subgraph subgraph) {
            Cost min_cost = subgraph.calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map());
            subgraphs.emplace_back(min_cost, std::move(subgraph));
            return subgraph_iterators::IterationControl::Continue;
        });

        // Sort subgraphs with decreasing costs.
        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

        Cost packing_cost = 0;
        VertexPairMap<bool> covered_by_packing(m_edit_state->graph().size());
        std::vector<Subgraph> packing;
        std::vector<VertexPair> min_cost_vertex_pairs;

        // Insert forbidden subgraphs with decreasing minimum edit cost into the bound.
        // Only insert a subgraph if it does not share an editable vertex pair with a subgraph already in the bound.
        for (const auto&[cost, subgraph] : subgraphs) {

            // Check if the subgraph is adjacent to one already used in the bound.
            bool touches_bound = false;
            for (VertexPair uv : subgraph.vertex_pairs())
                if (!m_edit_state->is_marked(uv) && covered_by_packing[uv]) {
                    touches_bound = true;
                    break;
                }

            if (!touches_bound) {
                VertexPair xy{static_cast<Vertex>(-2), static_cast<Vertex>(-1)};
                Cost min_cost = invalid_cost;
                for (VertexPair uv : subgraph.vertex_pairs()) {
                    if (!m_edit_state->is_marked(uv)) {
                        covered_by_packing[uv] = true;
                        if (m_edit_state->cost(uv) < min_cost) {
                            min_cost = m_edit_state->cost(uv);
                            xy = uv;
                        }
                    }
                }

                packing_cost += cost;
                packing.push_back(subgraph);
                min_cost_vertex_pairs.push_back(xy);
            }
        }

        return {packing_cost, covered_by_packing, packing, min_cost_vertex_pairs};
    }

    template class SortedGreedy<Options::FSG::C4P4>;
    template class SortedGreedy<Options::FSG::P3>;
}
