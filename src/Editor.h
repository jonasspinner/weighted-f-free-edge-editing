//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_EDITOR_H
#define CONCEPT_EDITOR_H

#include <iostream>
#include <utility>

#include "graph/Graph.h"
#include "graph/VertexPairMap.h"
#include "Instance.h"
#include "Configuration.h"
#include "Statistics.h"

#include "interfaces/SelectorI.h"
#include "interfaces/LowerBoundI.h"
#include "interfaces/ConsumerI.h"
#include "interfaces/FinderI.h"

#include "finder/finder_utils.h"
#include "lower_bound/lower_bound_utils.h"
#include "selector/selector_utils.h"

#include "consumer/SubgraphStats.h"


class Editor {
private:
    Instance m_instance;
    VertexPairMap<bool> m_marked;

    std::unique_ptr<LowerBoundI> m_lower_bound;
    std::unique_ptr<SelectorI> m_selector;
    std::unique_ptr<SubgraphStats> m_subgraph_stats;

    std::vector<ConsumerI *> m_consumers;
    std::shared_ptr<FinderI> m_finder;

    std::vector<VertexPair> edits;
    bool m_found_solution;

    // The idea was/is to mark all vertex pairs which editing cost is larger than the current available cost in advance
    // for each recursion call. This could potentially increase the quality of packing based lower bounds as more vertex
    // pairs are marked and therefore subgraphs which have these marked pairs are not considered adjacent.
    class VertexPairPreMarker {
        std::vector<VertexPair> m_vertex_pairs;
        std::vector<size_t> m_indices;
        std::vector<bool> m_needs_reset;

        const VertexPairMap<Cost> &m_costs;
        const std::vector<ConsumerI *> &m_consumers;

        VertexPairMap<bool> &m_marked;
    public:
        VertexPairPreMarker(const Graph &graph, const VertexPairMap<Cost> &costs,
                            const std::vector<ConsumerI *> &consumers, VertexPairMap<bool> &marked) :
                                m_indices({0}), m_costs(costs), m_consumers(consumers), m_marked(marked) {
            m_vertex_pairs.reserve(graph.size() * (graph.size() - 1) / 2);
            for (VertexPair uv : graph.vertexPairs())
                m_vertex_pairs.push_back(uv);

            std::sort(m_vertex_pairs.begin(), m_vertex_pairs.end(),
                      [&](VertexPair uv, VertexPair xy) { return m_costs[uv] > m_costs[xy]; });

            m_needs_reset.resize(m_vertex_pairs.size());
        }

        void mark(Cost k) {
            // All vertex pairs would be marked. Skip the work as only one selector call is necessary.
            if (k < m_costs[m_vertex_pairs.back()]) {
                m_indices.push_back(m_indices.back());
                return;
            };

            size_t index = m_indices.back();
            while (index < m_vertex_pairs.size() && m_costs[m_vertex_pairs[index]] > k) {
                VertexPair uv = m_vertex_pairs[index];
                if (!m_marked[uv]) {
                    for (auto &c : m_consumers) c->before_mark(uv);
                    m_marked[uv] = true;
                    for (auto &c : m_consumers) c->after_mark(uv);
                    m_needs_reset[index] = true;
                }
                ++index;
            }
            m_indices.push_back(index);
        }

        void unmark() {
            size_t end = m_indices.back();
            m_indices.pop_back();
            size_t begin = m_indices.back();
            for  (size_t index = end - 1; begin <= index && index < end; --index) {
                if (m_needs_reset[index]) {
                    auto uv = m_vertex_pairs[index];
                    assert(m_marked[uv]);

                    m_marked[uv] = false;
                    for (auto &c : m_consumers) c->after_unmark(uv);

                    m_needs_reset[index] = false;
                }
            }
        }
    } m_ordered_vertex_pairs;

    class PreMarkerGuard {
        VertexPairPreMarker &m_marker;
        bool m_marked = false;
    public:
        explicit PreMarkerGuard(VertexPairPreMarker& marker) : m_marker(marker) {}
        void mark(Cost k) { assert(!m_marked); m_marker.mark(k); m_marked = true; }
        ~PreMarkerGuard() { if (m_marked) m_marker.unmark(); }
    };

    Statistics m_stats;

    Configuration m_config;
public:
    explicit Editor(Instance instance, Configuration config) :
            m_instance(std::move(instance)), m_marked(m_instance.graph.size()), m_found_solution(false),
            m_ordered_vertex_pairs(m_instance.graph, m_instance.costs, m_consumers, m_marked), m_config(std::move(config)) {

        m_finder = Finder::make(m_config.forbidden_subgraphs);
        m_subgraph_stats = std::make_unique<SubgraphStats>(m_finder, m_instance, m_marked);
        m_consumers.emplace_back(m_subgraph_stats.get());

        m_selector = Selector::make(m_config.selector, m_finder, m_instance, m_marked, *m_subgraph_stats);
        m_lower_bound = LowerBound::make(m_config.lower_bound, m_finder, m_instance, m_marked, *m_subgraph_stats, m_config);

        m_consumers.emplace_back(m_lower_bound.get());
        m_consumers.emplace_back(m_selector.get());
    }

    [[nodiscard]] Cost initial_lower_bound() const {
        auto k = std::numeric_limits<Cost>::max();
        for (auto &c : m_consumers) c->initialize(k);
        return m_lower_bound->calculate_lower_bound(k);
    }

    /**
     * Initializes statistics and consumers. Calls recursive editing function.
     *
     * @param k The maximimum allowed editing cost.
     * @param result_cb A callback which is called when a result is found (vector<VertexPair> -> void).
     * @param prune_cb A callback which is called when a branch is pruned ((Cost, Cost) -> void).
     * @return Whether a solution was found.
     */
    //template<typename ResultCallback, typename PrunedCallback>
    bool edit(Cost k,
            const std::function<void(const std::vector<VertexPair>&)>& result_cb,
            const std::function<void(Cost, Cost)>& prune_cb = [](Cost, Cost) {},
            const std::function<bool(Cost)> &call_cb = [](Cost) { return false; }) {
        auto result_cb2 = [&](const std::vector<VertexPair> &e) { result_cb(e); return !m_config.find_all_solutions; };

        // init stats
        m_stats = Statistics(-k / 2, k, 100);
        m_found_solution = false;

        for (auto &c : m_consumers) c->initialize(k);
        edit_recursive(k, result_cb2, prune_cb, call_cb);
        return m_found_solution;
    }

    [[nodiscard]] const Statistics &stats() const {
        return m_stats;
    }

private:
    /**
     * Perform an edit step with remaining edit cost k.
     *
     * @param k The remaining editing cost.
     * @param result_cb A callback which is called when a result is found (vector<VertexPair> -> bool). If it returns true, the execution is stopped early.
     * @param prune_cb A callback which is called when a branch is pruned ((Cost, Cost) -> void).
     * @return Whether the execution was stopped early
     */
    //template<typename ResultCallback, typename PrunedCallback>
    bool edit_recursive(Cost k,
            const std::function<bool(const std::vector<VertexPair>&)> &result_cb,
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

            if (edit_recursive(k - costs[uv], result_cb, prune_cb, call_cb)) return_value = true;

            for (auto &c : m_consumers) c->pop_state();

            unedit_edge(uv);

            if (return_value) break;
        }

        // Iterating in reverse order because SubgraphStats keeps the history on a stack.
        for (auto uv = problem.pairs.rbegin(); uv != problem.pairs.rend(); ++uv)
            if (m_marked[*uv])
                unmark_edge(*uv);


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

#endif //CONCEPT_EDITOR_H
