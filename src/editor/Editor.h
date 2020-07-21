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

#include "VertexPairPreMarker.h"


class Editor {
private:
    Instance m_instance;
    VertexPairMap<bool> m_marked;
    VertexPairMap<bool> m_edited;

    std::unique_ptr<LowerBoundI> m_lower_bound;
    std::unique_ptr<SelectorI> m_selector;
    std::unique_ptr<ConsumerI> m_subgraph_stats;

    std::vector<ConsumerI *> m_consumers;

    std::vector<VertexPair> edits;
    bool m_found_solution;

    // The idea was/is to mark all vertex pairs which editing cost is larger than the current available cost in advance
    // for each recursion call. This could potentially increase the quality of packing based lower bounds as more vertex
    // pairs are marked and therefore subgraphs which have these marked pairs are not considered adjacent.
    VertexPairPreMarker m_ordered_vertex_pairs;

    Statistics m_stats;

    Configuration m_config;
public:
    explicit Editor(Instance instance, Configuration config) :
            m_instance(std::move(instance)), m_marked(m_instance.graph.size()), m_edited(m_instance.graph.size()),
            m_found_solution(false), m_ordered_vertex_pairs(m_instance.graph, m_instance.costs, m_consumers, m_marked),
            m_config(std::move(config)) {

        switch (config.forbidden_subgraphs) {
            case Options::FSG::C4P4:
                init_consumers<Options::FSG::C4P4>();
                break;
            default:
                throw std::runtime_error("Cannot initialize consumers for given set of forbidden subgraphs.");
        }
    }

    template<Options::FSG FSG>
    void init_consumers() {
        auto stats = std::make_unique<SubgraphStats<FSG>>(m_instance, m_marked);
        m_selector = selector::make<FSG>(m_config.selector, m_instance, m_marked, *stats);
        m_lower_bound = lower_bound::make<FSG>(m_config.lower_bound, m_instance, m_marked, *stats, m_config);

        m_subgraph_stats = std::move(stats);

        m_consumers.push_back(m_subgraph_stats.get());
        m_consumers.push_back(m_lower_bound.get());
        m_consumers.push_back(m_selector.get());
    }

    [[nodiscard]] Cost initial_lower_bound() const {
        auto k = std::numeric_limits<Cost>::max();
        for (auto &c : m_consumers) c->initialize(k);
        return m_lower_bound->calculate_lower_bound(k);
    }

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
        auto result_cb2 = [&](const std::vector<VertexPair> &e) {
            result_cb(e);
            return !m_config.find_all_solutions;
        };

        // init stats
        m_stats = Statistics(-k / 2, k, 100);
        m_found_solution = false;

        for (auto &c : m_consumers) c->initialize(k);
        switch (m_selector->recursion_type()) {
            case SelectorI::RecursionType::Subgraph:
                edit_recursive_subgraph(k, result_cb2, prune_cb, call_cb);
                break;
            case SelectorI::RecursionType::VertexPair:
                edit_recursive_vertex_pair(k, result_cb2, prune_cb, call_cb);
                break;
            default:
                throw std::runtime_error("Invalid recursion type");
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
    //template<typename ResultCallback, typename PrunedCallback>
    bool edit_recursive_subgraph(Cost k,
                        const std::function<bool(const std::vector<VertexPair> &)> &result_cb,
                        const std::function<void(Cost, Cost)> &prune_cb,
                        const std::function<bool(Cost)> &call_cb) {
        const VertexPairMap<Cost> &costs = m_instance.costs;

        m_stats.calls(k)++;
        if (call_cb(k)) return false;

        PreMarkerGuard guard(m_ordered_vertex_pairs);
        if (m_config.pre_mark_vertex_pairs)
            guard.mark(k);

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
            return result_cb(edits); // output graph
        }

        // recurse on problem pairs. keep vertex pairs marked between calls.
        bool return_value = false;
        for (VertexPair uv : problem.pairs) {
            assert(!m_marked[uv]);

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
            if (m_marked[*uv])
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
    //template<typename ResultCallback, typename PrunedCallback>
    bool edit_recursive_vertex_pair(Cost k,
                        const std::function<bool(const std::vector<VertexPair> &)> &result_cb,
                        const std::function<void(Cost, Cost)> &prune_cb,
                        const std::function<bool(Cost)> &call_cb) {
        const VertexPairMap<Cost> &costs = m_instance.costs;

        m_stats.calls(k)++;
        if (call_cb(k)) return false;

        PreMarkerGuard guard(m_ordered_vertex_pairs);
        if (m_config.pre_mark_vertex_pairs)
            guard.mark(k);

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
            return result_cb(edits); // output graph
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
        assert(!m_marked[uv]);
        for (auto &c : m_consumers) c->before_mark(uv);

        m_marked[uv] = true;

        for (auto &c : m_consumers) c->after_mark(uv);
    }

    /**
     * Edit the vertex pair. Call before_edit and after_edit.
     *
     * @param uv
     */
    void edit_edge(VertexPair uv) {
        assert(m_marked[uv]);
        Graph &G = m_instance.graph;

        for (auto &c : m_consumers) c->before_edit(uv);

        G.toggleEdge(uv);
        edits.push_back(uv);

        m_edited[uv] = true;

        for (auto &c : m_consumers) c->after_edit(uv);
    }

    /**
     * Undo the edit in the graph. Call before_unedit and after_unedit.
     *
     * @param uv
     */
    void unedit_edge(VertexPair uv) {
        assert(m_marked[uv]);
        Graph &G = m_instance.graph;

        for (auto &c : m_consumers) c->before_unedit(uv);

        G.toggleEdge(uv);
        edits.pop_back();

        m_edited[uv] = false;

        for (auto &c : m_consumers) c->after_unedit(uv);
    }

    /**
     * Unmark the vertex pair. Call after_unmark.
     *
     * @param uv
     */
    void unmark_edge(VertexPair uv) {
        assert(m_marked[uv]);
        m_marked[uv] = false;

        for (auto &c : m_consumers) c->after_unmark(uv);
    }
};

#endif //WEIGHTED_F_FREE_EDGE_EDITING_EDITOR_H
