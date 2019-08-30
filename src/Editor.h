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

#include "finder/NaiveP3.h"
#include "finder/NaiveC4P4.h"
#include "finder/CenterC4P4.h"
#include "finder/CenterP3.h"

#include "lower_bound/TrivialLowerBound.h"
#include "lower_bound/LocalSearchLowerBound.h"
#include "lower_bound/GreedyLowerBound.h"

#include "selector/FirstEditable.h"
#include "selector/LeastWeight.h"
#include "selector/MostMarkedPairs.h"

#include "consumer/SubgraphStats.h"


/**
 * ConsumerI::push and ConsumerI::pop operations bound on lifetime of Level object.
 */
class Level {
    const std::vector<ConsumerI *> &consumers;
public:
    Level(const std::vector<ConsumerI *> &consumers_, Cost k) : consumers(consumers_) {
        for (auto &consumer : consumers)
            consumer->push(k);
    }
    ~Level() {
        for (auto &consumer : consumers)
            consumer->pop();
    }
};


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

    Statistics m_stats;

public:
    explicit Editor(Instance instance, Options::Selector selector, Options::FSG forbidden, Options::LB lower_bound) :
            m_instance(std::move(instance)), m_marked(m_instance.graph.size()) {

        m_finder = make_finder(forbidden, m_instance);
        m_subgraph_stats = std::make_unique<SubgraphStats>(m_finder, m_instance, m_marked);
        // m_consumers.emplace_back(m_subgraph_stats.get());

        m_selector = make_selector(selector, m_finder, m_instance, m_marked);
        m_lower_bound = make_lower_bound(lower_bound, m_finder, m_instance, m_marked);

        m_consumers.emplace_back(m_lower_bound.get());
        m_consumers.emplace_back(m_selector.get());
    }

    template<typename ResultCallback, typename PrunedCallback>
    bool edit(Cost k, ResultCallback result, PrunedCallback pruned) {
        // init stats
        m_stats = Statistics(-k / 10, k, 10);
        return edit_r(k, result, pruned);
    }

    [[nodiscard]] const Statistics &stats() const {
        return m_stats;
    }

private:
    /**
     * Perform an edit step with remaining edit cost k.
     *
     * @param k
     * @return
     */
    template<typename ResultCallback, typename PrunedCallback>
    bool edit_r(Cost k, ResultCallback result, PrunedCallback pruned) {
        const VertexPairMap<Cost> &costs = m_instance.costs;

        m_stats.calls(k)++;

        Level L(m_consumers, k);

        /*std::cout << "marked:";
        for (VertexPair uv : m_instance.graph.vertexPairs())
            if (m_marked[uv]) std::cout << " " << uv;
        std::cout << "\n";*/


        //Cost sum = 0;
        //for (unsigned i = 0; i < edits.size(); ++i) std::cout << "  ";
        //for (VertexPair uv : edits) { std::cout << uv << " "; sum += costs[uv]; }
        //std::cout << sum << std::endl;

        auto lb = m_lower_bound->result(k);
        if (k < lb) {
            //for (unsigned i = 0; i < edits.size(); ++i) std::cout << "  ";
            pruned(k, lb);
            m_stats.prunes(k)++;
            return false; /* unsolvable, too few edits remaining */
        }

        auto problem = m_selector->result(k);
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

        std::vector<VertexPair> M;

        bool solved = false;
        for (VertexPair uv : problem.pairs) {
            if (m_marked[uv]) continue;

            // std::cout << "edit " << uv << "\n";
            mark_and_edit_edge(uv);
            M.push_back(uv);

            // std::cout << "edited " << uv << "\n";
            if (edit_r(k - costs[uv], result, pruned)) solved = true;

            unedit_edge(uv);
            // std::cout << "unedit " << uv << "\n";

            // if (solved) break;
        }

        for (VertexPair uv : M)
            unmark_edge(uv);

        return solved;
    }

    void mark_and_edit_edge(VertexPair uv) {
        assert(!m_marked[uv]);
        Graph &G = m_instance.graph;

        for (auto &c : m_consumers) c->before_mark_and_edit(uv);  // all
        for (auto &c : m_consumers) c->before_mark(uv);           // all

        m_marked[uv] = true;

        for (auto &c : m_consumers) c->after_mark(uv);            // subgraph_stats
        for (auto &c : m_consumers) c->before_edit(uv);           // subgraph_stats

        G.toggle_edge(uv);
        edits.push_back(uv);

        for (auto &c : m_consumers) c->after_edit(uv);            // subgraph_stats
        for (auto &c : m_consumers) c->after_mark_and_edit(uv);   // all
    }

    void unedit_edge(VertexPair uv) {
        assert(m_marked[uv]);
        Graph &G = m_instance.graph;

        for (auto &c : m_consumers) c->before_edit(uv);           // subgraph_stats

        G.toggle_edge(uv);
        edits.pop_back();

        for (auto &c : m_consumers) c->after_edit(uv);            // subgraph_stats
    }

    void unmark_edge(VertexPair uv) {
        assert(m_marked[uv]);
        m_marked[uv] = false;

        for (auto &c : m_consumers) c->after_unmark(uv);          // subgraph_stats
    }


    static std::shared_ptr<FinderI> make_finder(Options::FSG forbidden, const Instance &instance) {
        switch (forbidden) {
            case Options::FSG::P3:
                return std::make_shared<Finder::CenterP3>(instance.graph);
            case Options::FSG::P4C4:
                return std::make_shared<Finder::CenterC4P4>(instance.graph);
            default:
                assert(false);
                return nullptr;
        }
    }

    static std::unique_ptr<SelectorI>
    make_selector(Options::Selector selector, const std::shared_ptr<FinderI> &finder, const Instance &instance,
                  const VertexPairMap<bool> &marked) {
        switch (selector) {
            case Options::Selector::LeastWeight:
                return std::make_unique<Selector::LeastWeight>(instance.costs, finder, marked);
            case Options::Selector::FirstEditable:
                return std::make_unique<Selector::FirstEditable>(finder, marked);
            default:
                assert(false);
                return nullptr;
        }
    }

    static std::unique_ptr<LowerBoundI>
    make_lower_bound(Options::LB lower_bound, const std::shared_ptr<FinderI> &finder, const Instance &instance,
                     const VertexPairMap<bool> &marked) {
        switch (lower_bound) {
            case Options::LB::No:
                return std::make_unique<LowerBound::TrivialLowerBound>(finder);
            case Options::LB::LocalSearch:
                return std::make_unique<LocalSearchLowerBound>(instance, marked, finder);
            case Options::LB::Greedy:
                return std::make_unique<LowerBound::GreedyLowerBound>(instance, marked, finder);
            default:
                assert(false);
                return nullptr;
        }
    }
};

#endif //CONCEPT_EDITOR_H
