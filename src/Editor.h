//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_EDITOR_H
#define CONCEPT_EDITOR_H

#include <iostream>

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

    /*
    // The idea was/is to mark all vertex pairs which editing cost is larger than the current available cost in advance
    // for each recursion call. This could potentially increase the quality of packing based lower bounds as more vertex
    // pairs are marked and therefore subgraphs which have these marked pairs are not considered adjacent.
    class OrderedVertexPairs {
        std::vector<VertexPair> m_vertex_pairs;
        std::vector<std::vector<VertexPair>> m_stack;
        std::vector<size_t> m_indices;
    public:
        OrderedVertexPairs(const Graph &graph, const VertexPairMap<Cost> &costs) : m_indices({0}) {
            m_vertex_pairs.reserve(graph.size() * (graph.size() - 1) / 2);
            for (VertexPair uv : graph.vertexPairs())
                m_vertex_pairs.push_back(uv);

            std::sort(m_vertex_pairs.begin(), m_vertex_pairs.end(),
                      [&](VertexPair uv, VertexPair xy) { return costs[uv] > costs[xy]; });

            for (VertexPair uv : m_vertex_pairs) {
                std::cerr << costs[uv] << " ";
            }
            std::cerr << "\n";
        }

        void push(Cost k, const VertexPairMap<Cost> &costs, VertexPairMap<bool> &marked, std::vector<ConsumerI *> &consumers) {
            size_t index = m_indices.back();
            m_stack.emplace_back();
            while (costs[m_vertex_pairs[index]] > k) {
                VertexPair uv = m_vertex_pairs[index];
                if (!marked[uv]) {
                    for (auto &c : consumers) c->before_mark(uv);
                    marked[uv] = true;
                    for (auto &c : consumers) c->after_mark(uv);
                    m_stack.back().push_back(uv);
                }
                ++index;
            }
            m_indices.push_back(index);
        }

        void pop(VertexPairMap<bool> &marked, std::vector<ConsumerI *> &consumers) {
            m_indices.pop_back();
            for (auto it = m_stack.back().rbegin(); it != m_stack.back().rend(); ++it) {
                VertexPair uv = *it;
                assert(marked[uv]);
                marked[uv] = false;
                for (auto &c : consumers) c->after_unmark(uv);
            }
            m_stack.pop_back();
        }
    } m_ordered_vertex_pairs;
     */

    Statistics m_stats;

    bool m_find_all_solutions;

public:
    explicit Editor(Instance instance, Options::Selector selector, Options::FSG forbidden, Options::LB lower_bound) :
            m_instance(std::move(instance)), m_marked(m_instance.graph.size()), m_found_solution(false), /*m_ordered_vertex_pairs(m_instance.graph, m_instance.costs),*/
            m_find_all_solutions(true) {

        m_finder = Finder::make(forbidden, m_instance.graph);
        m_subgraph_stats = std::make_unique<SubgraphStats>(m_finder, m_instance, m_marked);
        m_consumers.emplace_back(m_subgraph_stats.get());

        m_selector = Selector::make(selector, m_finder, m_instance, m_marked, *m_subgraph_stats);
        m_lower_bound = LowerBound::make(lower_bound, m_finder, m_instance, m_marked, *m_subgraph_stats);

        m_consumers.emplace_back(m_lower_bound.get());
        m_consumers.emplace_back(m_selector.get());
    }

    Cost initialize(Cost k) {
        for (auto &c : m_consumers) c->initialize(k);
        return m_lower_bound->calculate_lower_bound(std::numeric_limits<Cost>::max());
    }

    /**
     * Initializes statistics and consumers. Calls recursive editing function.
     *
     * @param k The maximimum allowed editing cost.
     * @param result_cb A callback which is called when a result is found (vector<VertexPair> -> void).
     * @param prune_cb A callback which is called when a branch is pruned ((Cost, Cost) -> void).
     * @return Whether a solution was found.
     */
    template<typename ResultCallback, typename PrunedCallback>
    bool edit(Cost k, ResultCallback result_cb, PrunedCallback prune_cb) {
        // init stats
        m_stats = Statistics(-k / 2, k, 100);
        m_found_solution = false;

        edit_recursive(k, [&](const std::vector<VertexPair> &e) { result_cb(e); return !m_find_all_solutions; }, prune_cb);
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
    template<typename ResultCallback, typename PrunedCallback>
    bool edit_recursive(Cost k, ResultCallback result_cb, PrunedCallback prune_cb) {
        const VertexPairMap<Cost> &costs = m_instance.costs;

        m_stats.calls(k)++;

        // m_ordered_vertex_pairs.push(k, costs, m_marked, m_consumers);

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
        } else if (k == 0) {
            // unsolved, no edits remaining
            prune_cb(0, 0);
            m_stats.prunes(k)++;
            return false;
        }


        // recurse on problem pairs. keep vertex pairs marked between calls.
        bool return_value = false;
        for (VertexPair uv : problem.pairs) {
            assert(!m_marked[uv]);

            for (auto &c : m_consumers) c->push_state(k); // next_state(state)

            mark_and_edit_edge(uv);

            if (edit_recursive(k - costs[uv], result_cb, prune_cb)) return_value = true;

            for (auto &c : m_consumers) c->pop_state(); // destroy next_state

            unedit_edge(uv);

            if (return_value) break;
        }

        // Iterating in reverse order because SubgraphStats keeps the history on a stack.
        for (auto uv = problem.pairs.rbegin(); uv != problem.pairs.rend(); ++uv)
            if (m_marked[*uv])
                unmark_edge(*uv);

        // m_ordered_vertex_pairs.pop(m_marked, m_consumers);

        return return_value;
    }

    /**
     * Mark the vertex pair uv and make the edit in the graph.
     * The consumers are informed before and after the actions are performed.
     *
     * @param uv
     */
    void mark_and_edit_edge(VertexPair uv) {
        assert(!m_marked[uv]);
        Graph &G = m_instance.graph;

        for (auto &c : m_consumers) c->before_mark_and_edit(uv);  // all - next_state
        for (auto &c : m_consumers) c->before_mark(uv);           // all - state

        m_marked[uv] = true;

        for (auto &c : m_consumers) c->after_mark(uv);            // all - state
        for (auto &c : m_consumers) c->before_edit(uv);           // subgraph_stats

        G.toggleEdge(uv);
        edits.push_back(uv);

        for (auto &c : m_consumers) c->after_edit(uv);            // subgraph_stats
        for (auto &c : m_consumers) c->after_mark_and_edit(uv);   // all - next_state
    }

    /**
     * Undo the edit in the graph.
     *
     * @param uv
     */
    void unedit_edge(VertexPair uv) {
        assert(m_marked[uv]);
        Graph &G = m_instance.graph;

        for (auto &c : m_consumers) c->before_unedit(uv);           // subgraph_stats

        G.toggleEdge(uv);
        edits.pop_back();

        for (auto &c : m_consumers) c->after_unedit(uv);            // subgraph_stats
    }

    /**
     * Unmark the vertex pair.
     *
     * @param uv
     */
    void unmark_edge(VertexPair uv) {
        assert(m_marked[uv]);
        m_marked[uv] = false;

        for (auto &c : m_consumers) c->after_unmark(uv);          // subgraph_stats
    }
};

#endif //CONCEPT_EDITOR_H
