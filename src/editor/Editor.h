#ifndef WEIGHTED_F_FREE_EDGE_EDITING_EDITOR_H
#define WEIGHTED_F_FREE_EDGE_EDITING_EDITOR_H

#include <iostream>
#include <utility>

#include "../graph/Graph.h"
#include "../graph/VertexPairMap.h"
#include "../Instance.h"
#include "../Configuration.h"
#include "../Statistics.h"

#include "../selector/SelectorI.h"
#include "../lower_bound/LowerBoundI.h"
#include "../consumer/ConsumerI.h"

#include "../lower_bound/lower_bound_utils.h"
#include "../selector/selector_utils.h"

#include "../consumer/SubgraphStats.h"
#include "EditState.h"


class Editor {
private:
    std::unique_ptr<EditState> m_edit_state;

    std::unique_ptr<LowerBoundI> m_lower_bound;
    std::unique_ptr<SelectorI> m_selector;
    std::unique_ptr<ConsumerI> m_subgraph_stats;

    std::vector<ConsumerI *> m_consumers;

    bool m_found_solution = false;

    Cost m_initial_lower_bound = 0;
    bool m_is_initialized = false;

    Statistics m_stats;

    Configuration m_config;
    std::size_t m_num_calls_to_edit{0};
public:
    explicit Editor(Graph graph, VertexPairMap<Cost> costs, Configuration config) :
            m_edit_state(std::make_unique<EditState>(std::move(graph), std::move(costs))),
            m_config(std::move(config)) {

        switch (config.forbidden_subgraphs) {
            case Options::FSG::C4P4:
                make_consumers<Options::FSG::C4P4>();
                break;
            case Options::FSG::P3:
                make_consumers<Options::FSG::P3>();
                break;
            default:
                throw std::runtime_error("Cannot initialize consumers for given set of forbidden subgraphs.");
        }
    }

private:
    template<Options::FSG FSG>
    void make_consumers() {
        auto stats = std::make_unique<SubgraphStats<FSG>>(m_edit_state.get());
        m_selector = selector::make<FSG>(m_config.selector, m_edit_state.get(), stats.get());
        m_lower_bound = lower_bound::make<FSG>(m_config.lower_bound, m_edit_state.get(), stats.get(), m_config);

        m_subgraph_stats = std::move(stats);

        m_consumers.push_back(m_subgraph_stats.get());
        m_consumers.push_back(m_lower_bound.get());
        m_consumers.push_back(m_selector.get());
    }

    static bool allows_multiple_calls_to_edit(Options::FPTSearchStrategy strategy) {
        return strategy != Options::FPTSearchStrategy::Fixed;
    }

    void initalize_consumers(Cost k) {
        if (allows_multiple_calls_to_edit(m_config.search_strategy)) {
            k = std::numeric_limits<Cost>::max();
        }
        for (auto &c : m_consumers) {
            c->initialize(k);
        }
        m_is_initialized = true;
    }

public:

    [[nodiscard]] Cost initial_lower_bound() {
        if (!m_is_initialized) {
            auto k = std::numeric_limits<Cost>::max();
            initalize_consumers(k);
            m_initial_lower_bound = m_lower_bound->calculate_lower_bound(k);
        }
        return m_initial_lower_bound;
    }

private:
    enum class CallbackControl : bool {
        Continue = 0,
        Break = 1,
    };

    constexpr CallbackControl break_if(bool condition) noexcept {
        if (condition) {
            return CallbackControl::Break;
        } else {
            return CallbackControl::Continue;
        }
    }
public:

    /**
     * Initializes statistics and consumers. Calls recursive editing function.
     *
     * @param k The maximum allowed editing cost.
     * @param result_cb A callback which is called when a result is found (vector<VertexPair> -> void).
     * @param prune_cb A callback which is called when a branch is pruned ((Cost, Cost) -> void).
     * @return Whether a solution was found.
     */
    //template<typename ResultCallback, typename PrunedCallback>
    bool edit(Cost k,
              const std::function<void(const std::vector<VertexPair> &)> &result_cb,
              const std::function<void(Cost, Cost)> &prune_cb = [](Cost, Cost) {},
              const std::function<bool(Cost)> &call_cb = [](Cost) { return false; }) {

        ++m_num_calls_to_edit;
        if (!allows_multiple_calls_to_edit(m_config.search_strategy) && m_num_calls_to_edit >= 2) {
            throw std::runtime_error("This configuration only allows one call to edit(...)");
        }

        bool stop_after_first_solution = !m_config.find_all_solutions;

        auto result_cb_wrapper = [&](const std::vector<VertexPair> &e) {
            result_cb(e);
            return break_if(stop_after_first_solution);
        };

        auto call_cb_wrapper = [&](Cost cost) {
            return break_if(call_cb(cost));
        };

        // init stats
        m_stats = Statistics(-k / 2, k, 100);
        m_found_solution = false;

        if (!m_is_initialized) {
            initalize_consumers(k);
        }

        // If we expect multiple calls to edit(...), we have to save the initial state.
        if (allows_multiple_calls_to_edit(m_config.search_strategy)) {
            for (auto c : m_consumers) {
                c->push_state(k);
            }
        }


        switch (m_selector->recursion_type()) {
            case SelectorI::RecursionType::Subgraph:
                edit_recursive_subgraph(k, result_cb_wrapper, prune_cb, call_cb_wrapper);
                break;
            case SelectorI::RecursionType::VertexPair:
                edit_recursive_vertex_pair(k, result_cb_wrapper, prune_cb, call_cb_wrapper);
                break;
            default:
                throw std::runtime_error("Invalid recursion type");
        }


        // Restore the initial state
        if (allows_multiple_calls_to_edit(m_config.search_strategy)) {
            for (auto c : m_consumers) {
                c->pop_state();
            }
        }

        return m_found_solution;
    }

    [[nodiscard]] const Statistics &stats() const {
        return m_stats;
    }

private:
    /**
     * Perform an edit step with remaining edit cost k. Branch on all unmarked vertex pairs of a single subgraphs.
     *
     * @param k The remaining editing cost.
     * @param result_cb A callback which is called when a result is found (vector<VertexPair> -> bool). If it returns true, the execution is stopped early.
     * @param prune_cb A callback which is called when a branch is pruned ((Cost, Cost) -> void).
     * @return Whether the execution was stopped early
     */
    template<typename ResultCallback, typename PrunedCallback, typename CallCallback>
    bool edit_recursive_subgraph(Cost k,
                                 const ResultCallback &result_cb,
                                 const PrunedCallback &prune_cb,
                                 const CallCallback &call_cb) {
        static_assert(std::is_invocable_r_v<CallbackControl, ResultCallback, const std::vector<VertexPair> &>,
                      "ResultCallback must have CallbackControl(const std::vector<VertexPair> &) signature.");
        static_assert(std::is_invocable_r_v<void, PrunedCallback, Cost, Cost>,
                      "PrunedCallback must have void(Cost, Cost) signature.");
        static_assert(std::is_invocable_r_v<CallbackControl, CallCallback, Cost>,
                      "CallCallback must have CallbackControl(Cost) signature.");

        const auto &costs = m_edit_state->cost_map();
        const auto &edits = m_edit_state->edits();

        m_stats.calls(k)++;
        if (call_cb(k) == CallbackControl::Break) return false;

        auto lb = m_lower_bound->calculate_lower_bound(k);
        if (k < lb) {
            // unsolvable, too few edits remaining
            prune_cb(k, lb);
            m_stats.prunes(k)++;
            return false;
        }

        auto problem = m_selector->select_problem(k);

        if (problem.solved) {
            // solved
            m_found_solution = true;
            return result_cb(edits) == CallbackControl::Break; // output graph
        }

        // recurse on problem pairs. keep vertex pairs marked between calls.
        bool return_value = false;
        for (VertexPair uv : problem.pairs) {
            assert(!m_edit_state->is_marked(uv));

            mark_edge(uv);

            for (auto &c : m_consumers) c->push_state(k);

            edit_edge(uv);

            if (edit_recursive_subgraph(k - costs[uv], result_cb, prune_cb, call_cb)) return_value = true;

            for (auto &c : m_consumers) c->pop_state();

            unedit_edge(uv);

            if (return_value) break;

            if (k < m_lower_bound->calculate_lower_bound_no_edit_branch()) break;
        }

        // Iterating in reverse order because SubgraphStats keeps the history on a stack.
        for (auto uv = problem.pairs.rbegin(); uv != problem.pairs.rend(); ++uv)
            if (m_edit_state->is_marked(*uv))
                unmark_edge(*uv);


        return return_value;
    }


    /**
     * Perform an edit step with remaining edit cost k. Branch on vertex pairs.
     *
     * @param k The remaining editing cost.
     * @param result_cb A callback which is called when a result is found (vector<VertexPair> -> bool). If it returns true, the execution is stopped early.
     * @param prune_cb A callback which is called when a branch is pruned ((Cost, Cost) -> void).
     * @return Whether the execution was stopped early
     */
    template<typename ResultCallback, typename PrunedCallback, typename CallCallback>
    bool edit_recursive_vertex_pair(Cost k,
                                    const ResultCallback &result_cb,
                                    const PrunedCallback &prune_cb,
                                    const CallCallback &call_cb) {
        static_assert(std::is_invocable_r_v<CallbackControl, ResultCallback, const std::vector<VertexPair> &>,
                      "ResultCallback must have CallbackControl(const std::vector<VertexPair> &) signature.");
        static_assert(std::is_invocable_r_v<void, PrunedCallback, Cost, Cost>,
                      "PrunedCallback must have void(Cost, Cost) signature.");
        static_assert(std::is_invocable_r_v<CallbackControl, CallCallback, Cost>,
                      "CallCallback must have CallbackControl(Cost) signature.");

        const auto &costs = m_edit_state->cost_map();
        const auto &edits = m_edit_state->edits();

        m_stats.calls(k)++;
        if (call_cb(k) == CallbackControl::Break) return false;

        auto lb = m_lower_bound->calculate_lower_bound(k);
        if (k < lb) {
            // unsolvable, too few edits remaining
            prune_cb(k, lb);
            m_stats.prunes(k)++;
            return false;
        }

        auto problem = m_selector->select_problem(k);

        if (problem.solved) {
            // solved
            m_found_solution = true;
            return result_cb(edits) == CallbackControl::Break; // output graph
        }


        bool return_value = false;

        if (!problem.pairs.empty()) {
            assert(problem.pairs.size() == 1);
            VertexPair uv = problem.pairs[0];


            mark_edge(uv);

            // uv not edited
            for (auto &c : m_consumers) c->push_state(k);

            if (edit_recursive_vertex_pair(k, result_cb, prune_cb, call_cb)) return_value = true;

            for (auto &c : m_consumers) c->pop_state();


            if (!return_value) {
                // uv edited
                for (auto &c : m_consumers) c->push_state(k);

                edit_edge(uv);

                if (edit_recursive_vertex_pair(k - costs[uv], result_cb, prune_cb, call_cb)) return_value = true;

                for (auto &c : m_consumers) c->pop_state();

                unedit_edge(uv);
            }


            unmark_edge(uv);
        }

        return return_value;
    }

    /**
     * Mark the vertex pair. Call before_mark and after_mark.
     *
     * @param uv
     */
    void mark_edge(VertexPair uv) {
        assert(!m_edit_state->is_marked(uv));
        for (auto &c : m_consumers) c->before_mark(uv);

        m_edit_state->mark_edge(uv);

        for (auto &c : m_consumers) c->after_mark(uv);
    }

    /**
     * Edit the vertex pair. Call before_edit and after_edit.
     *
     * @param uv
     */
    void edit_edge(VertexPair uv) {
        assert(m_edit_state->is_marked(uv));

        for (auto &c : m_consumers) c->before_edit(uv);

        m_edit_state->edit_edge(uv);

        for (auto &c : m_consumers) c->after_edit(uv);
    }

    /**
     * Undo the edit in the graph. Call before_unedit and after_unedit.
     *
     * @param uv
     */
    void unedit_edge(VertexPair uv) {
        assert(m_edit_state->is_marked(uv));

        for (auto &c : m_consumers) c->before_unedit(uv);

        m_edit_state->unedit_edge(uv);

        for (auto &c : m_consumers) c->after_unedit(uv);
    }

    /**
     * Unmark the vertex pair. Call after_unmark.
     *
     * @param uv
     */
    void unmark_edge(VertexPair uv) {
        assert(m_edit_state->is_marked(uv));

        m_edit_state->unmark_edge(uv);

        for (auto &c : m_consumers) c->after_unmark(uv);
    }
};

#endif //WEIGHTED_F_FREE_EDGE_EDITING_EDITOR_H
