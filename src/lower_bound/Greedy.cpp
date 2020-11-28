#include "Greedy.h"


namespace lower_bound {
    /**
     * Calculates a lower bound on the costs required to solve the current instance.
     *
     * Inserts forbidden subgraphs into the bound if possible. A subgraph is not inserted if it shares an editable
     * vertex pair with a subgraph already in the bound.
     *
     * Complexity:
     *      m = #forbidden subgraphs
     *      O(m) time
     *      O(1) additional space
     *
     * @param k Not used
     * @return A lower bound on the costs required to solve the current instance.
     */
    template <Options::FSG SetOfForbiddenSubgraphs>
    Cost Greedy<SetOfForbiddenSubgraphs>::calculate_lower_bound(Cost k) {
        Cost bound_size = 0;
        m_used_in_bound.clear();

        finder.find_unique(m_edit_state->graph(), [&](Subgraph subgraph) {

            // Check if the subgraph is adjacent to one already used in the bound.
            bool touches_bound = false;
            for (auto uv : subgraph.non_converting_edits()) {
                if (!m_edit_state->is_marked(uv) && m_used_in_bound[uv]) {
                    touches_bound = true;
                    break;
                }
            }

            if (!touches_bound) {
                Cost cost = subgraph.calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map());
                if (cost == invalid_cost) {
                    bound_size = invalid_cost;
                    return subgraph_iterators::IterationControl::Break;
                }
                bound_size += cost;

                for (auto uv : subgraph.non_converting_edits()) {
                    if (!m_edit_state->is_marked(uv))
                        m_used_in_bound[uv] = true;
                }
            }
            return subgraph_iterators::break_if(bound_size > k);
        });

        return bound_size;
    }

    template class Greedy<Options::FSG::C4P4>;
}
