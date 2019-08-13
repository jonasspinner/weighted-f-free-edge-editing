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

#include "interfaces/SelectorI.h"
#include "interfaces/LowerBoundI.h"
#include "interfaces/ConsumerI.h"
#include "interfaces/FinderI.h"

#include "finder/NaiveP3.h"
#include "finder/NaiveC4P4.h"
#include "finder/CenterC4P4.h"
#include "finder/CenterP3.h"

#include "lower_bound/NoLowerBound.h"
#include "lower_bound/IteratedLocalSearch.h"
#include "lower_bound/GreedyLowerBound.h"

#include "selector/FirstEditable.h"
#include "selector/LeastWeight.h"

#include "consumer/SubgraphStats.h"


#define call_for_each(c, m, s, e) \
do { for (size_t i = 0; i < c.size(); ++i) {\
    c[i]->m(*s[i], e); \
}} while(0)


class Editor {
private:
    using States = std::vector<std::unique_ptr<StateI>>;

    Instance m_instance;
    VertexPairMap<bool> m_marked;
    std::unique_ptr<LowerBoundI> m_lower_bound;
    std::unique_ptr<SelectorI> m_selector;
    std::vector<ConsumerI *> m_consumers;
    std::shared_ptr<FinderI> m_finder;
    std::unique_ptr<SubgraphStats> m_subgraph_stats;

    std::vector<VertexPair> edits;

public:
    explicit Editor(Instance instance, Options::Selector selector, Options::FSG forbidden, Options::LB lower_bound) :
            m_instance(std::move(instance)), m_marked(m_instance.graph.size()) {

        m_finder = make_finder(forbidden, m_instance);
        m_selector = make_selector(selector, m_finder, m_instance, m_marked);
        m_lower_bound = make_lower_bound(lower_bound, m_finder, m_instance, m_marked);

        m_subgraph_stats = std::make_unique<SubgraphStats>(m_finder, m_instance, m_marked);

        m_consumers.emplace_back(m_lower_bound.get());
        m_consumers.emplace_back(m_selector.get());

        // m_consumers.emplace_back(m_subgraph_stats.get());
    }

    template<typename ResultCallback, typename PrunedCallback>
    bool edit(Cost k, ResultCallback result, PrunedCallback pruned) {
        // init stats
        auto states = make_states(k);
        return edit_r(k, std::move(states), result, pruned);
    }

private:
    /**
     * Perform an edit step with remaining edit cost k.
     *
     * @param k
     * @return
     */
    template<typename ResultCallback, typename PrunedCallback>
    bool edit_r(Cost k, States states, ResultCallback result, PrunedCallback pruned) {
        const VertexPairMap<Cost> &costs = m_instance.costs;

        /*std::cout << "marked:";
        for (VertexPair uv : m_instance.graph.vertexPairs())
            if (m_marked[uv]) std::cout << " " << uv;
        std::cout << "\n";*/


        //Cost sum = 0;
        //for (unsigned i = 0; i < edits.size(); ++i) std::cout << "  ";
        //for (VertexPair uv : edits) { std::cout << uv << " "; sum += costs[uv]; }
        //std::cout << sum << std::endl;

        auto lb = m_lower_bound->result(*states[0], k);
        if (k < lb) {
            //for (unsigned i = 0; i < edits.size(); ++i) std::cout << "  ";
            pruned(k, lb);
            return false; /* unsolvable, too few edits remaining */
        }

        auto problem = m_selector->result(*states[1], k);
        if (problem.solved) {
            result(edits); // output graph
            return true; /* solved */
        } else if (k == 0) {
            pruned(0, 0);
            return false;
        } /* unsolved, no edits remaining */


        /* std::cout << "recurse on problem set {";
        for (auto uv : problem.pairs) std::cout << " " << uv;
        std::cout << " }\n"; */


        bool solved = false;
        for (VertexPair uv : problem.pairs) {
            if (m_marked[uv]) continue;

            // std::cout << "edit " << uv << "\n";
            mark_and_edit_edge(states, uv);

            auto next_states = copy(states);

            // std::cout << "edited " << uv << "\n";
            if (edit_r(k - costs[uv], std::move(next_states), result, pruned)) solved = true;

            unedit_edge(states, uv);
            // std::cout << "unedit " << uv << "\n";

            // if (solved) break;
        }

        for (VertexPair uv : problem.pairs) {
            if (m_marked[uv]) unmark_edge(states, uv);
        }

        return solved;
    }

    void mark_and_edit_edge(States &states, VertexPair uv) {
        assert(!m_marked[uv]);
        Graph &G = m_instance.graph;

        call_for_each(m_consumers, before_mark_and_edit, states, uv);   // all
        call_for_each(m_consumers, before_mark, states, uv);            // all

        m_marked[uv] = true;

        call_for_each(m_consumers, after_mark, states, uv);             // subgraph_stats
        call_for_each(m_consumers, before_edit, states, uv);            // subgraph_stats

        G.toggle_edge(uv);
        edits.push_back(uv);

        call_for_each(m_consumers, after_edit, states, uv);             // subgraph_stats
        call_for_each(m_consumers, after_mark_and_edit, states, uv);    // all

    }

    void unedit_edge(States &states, VertexPair uv) {
        assert(m_marked[uv]);
        Graph &G = m_instance.graph;

        call_for_each(m_consumers, before_edit, states, uv);            // subgraph_stats

        G.toggle_edge(uv);
        edits.pop_back();

        call_for_each(m_consumers, after_edit, states, uv);             // subgraph_stats

    }

    void unmark_edge(States &states, VertexPair uv) {
        assert(m_marked[uv]);
        m_marked[uv] = false;

        call_for_each(m_consumers, after_unmark, states, uv);           // subgraph_stats
    }

    States make_states(Cost k) {
        States states;
        for (auto consumer : m_consumers) {
            states.push_back(consumer->initialize(k));
        }
        return states;
    }

    static States copy(const States &states) {
        States new_states;
        for (const auto &state : states) {
            new_states.push_back(state->copy());
        }
        return new_states;
    }

    static std::shared_ptr<FinderI> make_finder(Options::FSG forbidden, const Instance &instance) {
        switch (forbidden) {
            case Options::FSG::P3:
                return std::make_shared<Finder::CenterP3>(instance.graph);
            case Options::FSG::P4C4:
                return std::make_shared<Finder::CenterC4P4>(instance.graph);
        }
        assert(false);
    }

    static std::unique_ptr<SelectorI>
    make_selector(Options::Selector selector, const std::shared_ptr<FinderI> &finder, const Instance &instance,
                  const VertexPairMap<bool> &marked) {
        switch (selector) {
            case Options::Selector::LeastWeight:
                return std::make_unique<Selector::LeastWeight>(instance.costs, finder, marked);
            case Options::Selector::FirstEditable:
                return std::make_unique<Selector::FirstEditable>(finder, marked);
        }
        assert(false);
    }

    static std::unique_ptr<LowerBoundI>
    make_lower_bound(Options::LB lower_bound, const std::shared_ptr<FinderI> &finder, const Instance &instance,
                     const VertexPairMap<bool> &marked) {
        switch (lower_bound) {
            case Options::LB::No:
                return std::make_unique<LowerBound::NoLowerBound>(finder);
            case Options::LB::LocalSearch:
                return std::make_unique<IteratedLocalSearch>(instance, marked, finder);
            case Options::LB::Greedy:
                return std::make_unique<LowerBound::GreedyLowerBound>(instance, marked, finder);
        }
        assert(false);
    }
};

#undef call_for_each

#endif //CONCEPT_EDITOR_H
