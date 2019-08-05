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
#include "finder/CenterC4P4.h"
#include "finder/CenterP3.h"

#include "lower_bound/NoLowerBound.h"

#include "selector/FirstEditable.h"
#include "selector/LeastWeight.h"
#include "lower_bound/IteratedLocalSearch.h"
#include "lower_bound/GreedyLowerBound.h"
#include "consumer/SubgraphStats.h"


class Editor {
private:
    using LBI = LowerBoundI;
    using LB = LowerBound::NoLowerBound;
    using SelI = SelectorI;
    using Sel = Selector::LeastWeight;
    using CI = ConsumerI;

    using States = std::vector<std::unique_ptr<StateI>>;

    Instance m_instance;
    VertexPairMap<bool> m_forbidden;
    std::unique_ptr<LowerBoundI> m_lower_bound;
    std::unique_ptr<SelectorI> m_selector;
    std::vector<ConsumerI*> m_consumers;
    std::shared_ptr<FinderI> m_finder;
    std::unique_ptr<SubgraphStats> m_subgraph_stats;

    std::vector<VertexPair> edits;

public:
    explicit Editor(Instance instance, Configuration::SelectorOption selector, Configuration::ForbiddenSubgraphs forbidden, Configuration::LowerBound lower_bound) : m_instance(std::move(instance)), m_forbidden(m_instance.graph) {
        switch (forbidden) {
            case Configuration::ForbiddenSubgraphs::P3:
                m_finder = std::make_shared<Finder::CenterP3>(m_instance.graph);
                break;
            case Configuration::ForbiddenSubgraphs::P4C4:
                m_finder = std::make_shared<Finder::CenterC4P4>(m_instance.graph);
                break;
        }
        switch (selector) {
            case Configuration::SelectorOption::LeastWeight:
                m_selector = std::make_unique<Selector::LeastWeight>(m_instance.costs, m_finder, m_forbidden);
                break;
            case Configuration::SelectorOption::FirstEditable:
                m_selector = std::make_unique<Selector::FirstEditable>(m_instance.graph, m_instance.costs, m_finder, m_forbidden);
                break;
        }
        switch (lower_bound) {
            case Configuration::LowerBound::No:
                m_lower_bound = std::make_unique<LowerBound::NoLowerBound>(m_finder);
                break;
            case Configuration::LowerBound::IteratedLocalSearch:
                m_lower_bound = std::make_unique<IteratedLocalSearch>(m_instance, m_forbidden, m_finder);
                break;
            case Configuration::LowerBound::Greedy:
                m_lower_bound = std::make_unique<LowerBound::GreedyLowerBound>(m_instance, m_forbidden, m_finder);
                break;
        }
        m_lower_bound = std::make_unique<LowerBound::GreedyLowerBound>(m_instance, m_forbidden, m_finder);
        m_subgraph_stats = std::make_unique<SubgraphStats>(m_finder, m_instance, m_forbidden);

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
        // std::cout << "edit(" << k << ")\n";
        auto lb = m_lower_bound->result(*states[0], k);
        if (k < lb) {
            // std::cout << "pruned\tlb(" << k << ") = " << lb << ", eps = " << lb - k << "\n";
            pruned(k, lb);
            return false; /* unsolvable, too few edits remaining */
        }

        auto problem = m_selector->result(*states[1], k);
        if (problem.solved) {
            // output graph
            result(edits);
            return true; /* solved */
        } else if (k == 0) { pruned(0, 0); return false; } /* unsolved, no edits remaining */


        /* std::cout << "recurse on problem set {";
        for (auto uv : problem.pairs) std::cout << " " << uv;
        std::cout << " }\n"; */


        bool solved = false;
        for (auto uv : problem.pairs) {
            if (m_forbidden[uv]) continue;

            auto next_states = copy(states);

            // std::cout << "edit " << uv << "\n";
            mark_and_edit_edge(states, uv);

            // std::cout << "edited " << uv << "\n";
            if (edit_r(k - costs[uv], std::move(next_states), result, pruned)) solved = true;

            unedit_edge(states, uv);
            // std::cout << "unedit " << uv << "\n";

            // if (solved) break;
        }

        for (const auto &uv : problem.pairs) {
            if (m_forbidden[uv]) unmark_edge(states, uv);
        }

        return solved;
    }

    void mark_and_edit_edge(States &states, const VertexPair &uv) {
        assert(!m_forbidden[uv]);
        Graph &G = m_instance.graph;
        // before_mark_and_edit         all
        for (size_t i = 0; i < m_consumers.size(); ++i) {
            m_consumers[i]->before_mark_and_edit(*states[i], uv);
        }
        // before_mark                  all
        for (size_t i = 0; i < m_consumers.size(); ++i) {
            m_consumers[i]->before_mark(*states[i], uv);
        }

        m_forbidden[uv] = true;

        // after_mark                   subgraph_stats
        for (size_t i = 0; i < m_consumers.size(); ++i) {
            m_consumers[i]->after_mark(*states[i], uv);
        }
        // before_edit                  subgraph_stats
        for (size_t i = 0; i < m_consumers.size(); ++i) {
            m_consumers[i]->before_edit(*states[i], uv);
        }

        G.toggle_edge(uv);
        edits.push_back(uv);

        // after_edit                   subgraph_stats
        for (size_t i = 0; i < m_consumers.size(); ++i) {
            m_consumers[i]->after_edit(*states[i], uv);
        }
        // after_mark_and_edit          all
        for (size_t i = 0; i < m_consumers.size(); ++i) {
            m_consumers[i]->after_mark_and_edit(*states[i], uv);
        }
    }

    void unedit_edge(States &states, const VertexPair &uv) {
        assert(m_forbidden[uv]);
        Graph &G = m_instance.graph;

        // before_edit                  subgraph_stats
        for (size_t i = 0; i < m_consumers.size(); ++i) {
            m_consumers[i]->before_edit(*states[i], uv);
        }

        G.toggle_edge(uv);
        edits.pop_back();

        // after_edit                   subgraph_stats
        for (size_t i = 0; i < m_consumers.size(); ++i) {
            m_consumers[i]->after_edit(*states[i], uv);
        }
    }

    void unmark_edge(States &states, const VertexPair &uv) {
        assert(m_forbidden[uv]);
        m_forbidden[uv] = false;

        // after_unmark                 subgraph_stats
        for (size_t i = 0; i < m_consumers.size(); ++i) {
            m_consumers[i]->after_unmark(*states[i], uv);
        }
    }

    std::vector<std::unique_ptr<StateI>> make_states(Cost k) {
        std::vector<std::unique_ptr<StateI>> states;
        for (auto consumer : m_consumers) {
            states.push_back(consumer->initialize(k));
        }
        return states;
    }

    static States copy(const States& states) {
        States new_states;
        for (const auto& state : states) {
            new_states.push_back(state->copy());
        }
        return new_states;
    }
};


#endif //CONCEPT_EDITOR_H
