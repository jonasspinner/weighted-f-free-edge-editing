//
// Created by jonas on 25.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LOCALSEARCHLOWERBOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LOCALSEARCHLOWERBOUND_H


#include <random>
#include "../interfaces/LowerBoundI.h"
#include "../interfaces/FinderI.h"

class LocalSearchLowerBound : public LowerBoundI {
private:
    class State {
    public:
        struct Element {
            Cost cost;
            Subgraph subgraph;
        };
        std::vector<Element> bound;
        Cost m_cost = 0;

        [[nodiscard]] Cost cost() const { return m_cost; }

        void remove(size_t index) {
            m_cost -= bound[index].cost;
            bound[index] = std::move(bound.back());
            bound.pop_back();
        }

        void replace(size_t index, Element &&element) {
            m_cost -= bound[index].cost;
            m_cost += element.cost;
            bound[index] = std::move(element);
        }

        void insert(Element &&element) {
            m_cost += element.cost;
            bound.emplace_back(std::move(element));
        }

        void set_unsolvable() {
            bound.clear();
            m_cost = invalid_cost;
        }

        [[nodiscard]] bool solvable() const { return m_cost != invalid_cost; }

        void clear() {
            bound.clear();
            m_cost = 0;
        }

        void recalculate(const VertexPairMap<bool> &marked, const VertexPairMap<Cost> &costs) {
            m_cost = 0;
            for (auto &e : bound) {
                e.cost = ::cost(e.subgraph, marked, costs);
                m_cost += e.cost;
                if (e.cost == invalid_cost) {
                    set_unsolvable();
                    break;
                }
            }

#ifndef NDEBUG
            if (solvable()) {
                Cost sum = 0;
                for (const auto &e : bound) {
                    assert(e.cost == ::cost(e.subgraph, marked, costs));
                    sum += e.cost;
                }
                assert(m_cost == sum);
            }
#endif
        }

        void initialize_bound_graph(const VertexPairMap<bool> &marked, Graph &bound_graph) {
            bound_graph.clear_edges();

            for (const auto &[cost, subgraph] : bound) {
                for (VertexPair uv : subgraph.vertexPairs()) {
                    if (!marked[uv]) {
                        assert(!bound_graph.has_edge(uv));
                        bound_graph.set_edge(uv);
                    }
                }
            }

#ifndef NDEBUG
            VertexPairMap<bool> debug(bound_graph.size());
            for (const auto&[cost, subgraph] : bound)
                for (VertexPair uv : subgraph.vertexPairs())
                    if (!marked[uv]) {
                        assert(bound_graph.has_edge(uv));
                        assert(!debug[uv]);
                        debug[uv] = true;
                    }

            for (VertexPair uv : bound_graph.vertexPairs())
                if (marked[uv]) {
                    assert(!bound_graph.has_edge(uv));
                    assert(!debug[uv]);
                } else {
                    assert(bound_graph.has_edge(uv) == static_cast<bool>(debug[uv]));
                }
#endif
        }

        void assert_valid(const VertexPairMap<bool> &marked, const Graph &bound_graph, const VertexPairMap<Cost> &costs) {
            if (!solvable()) return;
            // every subgraph is valid
            //for (const auto &[_, subgraph] : bound) {
            //    verify(subgraph, graph);
            //}

            // costs are matching
            Cost sum = 0;
            for (const auto &e : bound) {
                assert(e.cost == ::cost(e.subgraph, marked, costs));
                sum += e.cost;
            }
            assert(m_cost == sum);

            VertexPairMap<bool> debug(bound_graph.size());

            for (const auto& [cost, subgraph] : bound) {
                for (VertexPair uv : subgraph.vertexPairs()) {
                    if (!marked[uv]) {
                        assert(bound_graph.has_edge(uv));
                        assert(!debug[uv]);
                        debug[uv] = true;
                    }
                }
            }

            for (VertexPair uv : bound_graph.vertexPairs()) {
                if (marked[uv]) {
                    assert(!bound_graph.has_edge(uv));
                    assert(!debug[uv]);
                } else {
                    assert(bound_graph.has_edge(uv) == static_cast<bool>(debug[uv]));
                }
            }

        }

        friend std::ostream &operator<<(std::ostream &os, const State &state) {
            if (!state.solvable()) return os << "unsolvable";
            std::cout << state.cost() << ":";
            for (const auto &[cost, subgraph] : state.bound)
                std::cout << " (" << subgraph << ", " << cost << ")";
            return os;
        }
    };

    Graph m_bound_graph;
    const VertexPairMap<Cost> &m_costs;
    const VertexPairMap<bool> &m_marked;
    std::vector<std::unique_ptr<State>> states;

public:
    explicit LocalSearchLowerBound(const Instance &instance,
                                   const VertexPairMap<bool> &forbidden, std::shared_ptr<FinderI> finder_ref) :
            LowerBoundI(std::move(finder_ref)), m_bound_graph(instance.graph.size()),
            m_costs(instance.costs), m_marked(forbidden) {}

    /**
     *
     * @param k
     * @return A lower bound on the costs required to solve the current instance.
     */
    Cost result(Cost k) override {
        auto &state = current_state();
        if (!state.solvable()) return state.cost();

        // std::cout << "result start: " << state << "\n";

        state.recalculate(m_marked, m_costs);
        state.initialize_bound_graph(m_marked, m_bound_graph);
        state.assert_valid(m_marked, m_bound_graph, m_costs);

        if (state.cost() <= k) {
            optimize_bound(state, k);
        }

        assert(state.cost() >= 0);

        // std::cout << "result end: " << state << "\n";
        return state.cost();
    }

    /**
     * Initializes the state by greedily constructing a maximal lower bound.
     *
     * @return
     */
    void initialize() override {
        using Element = State::Element;
        states.push_back(std::make_unique<State>());
        State &state = *states.back();

        m_bound_graph.clear_edges();

        std::vector<Element> subgraphs;

        finder->find([&](Subgraph &&subgraph) {
            Cost min_cost = cost(subgraph, m_marked, m_costs);
            subgraphs.push_back({min_cost, std::move(subgraph)});
            return false;
        });

        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const Element &lhs, const Element &rhs) { return lhs.cost > rhs.cost; });

        if (!subgraphs.empty() && subgraphs[0].cost == invalid_cost) {
            state.set_unsolvable();
            return;
        }

        for (auto& [cost, subgraph] : subgraphs) {
            bool inserted = try_insert_into_graph(subgraph, m_marked, m_bound_graph);

            if (inserted) {
                assert(cost != invalid_cost);
                state.insert({cost, std::move(subgraph)});
            }
        }
    }

    /**
     * Copy the state from the previous level.
     *
     * @param k
     * @return
     */
    void push_state(Cost /*k*/) override {
        assert(!states.empty());
        states.push_back(std::make_unique<State>(*states.back()));
    }

    void pop_state() override {
        assert(!states.empty());
        states.pop_back();
    }

    State &current_state() {
        assert(!states.empty());
        return *states.back();
    }

    State &parent_state() {
        assert(states.size() >= 2);
        return *states[states.size() - 2];
    }

    /**
     * Removes subgraphs from state.bound which have the vertex pair uv.
     *
     * @param uv
     */
    void before_mark_and_edit(VertexPair uv) override {
        auto &state = current_state();
        assert(state.solvable());
        // state.assert_valid(m_marked, m_bound_graph, m_costs);

        for (size_t i = 0; i < state.bound.size(); ) {
            const auto &[cost, subgraph] = state.bound[i];
            const auto vertices = subgraph.vertices();

            bool has_u = std::any_of(vertices.begin(), vertices.end(), [&](Vertex x) { return x == uv.u; });
            bool has_v = std::any_of(vertices.begin(), vertices.end(), [&](Vertex x) { return x == uv.v; });

            if (has_u && has_v) {
                // std::cout << "before_mark_and_edit remove " << subgraph;
                state.remove(i);
                // std::cout << "\t => " << state << "\n";
            } else {
                ++i;
            }
        }
    }

    void after_mark_and_edit(VertexPair uv) override {
        auto &state = current_state();
        assert(state.solvable());

        // TODO: Check if needed
        state.initialize_bound_graph(m_marked, m_bound_graph);

        std::vector<std::pair<Cost, Subgraph>> subgraphs;

        // The finder iterates over subgraphs having u and v as vertices.
        finder->find_near(uv, m_bound_graph, [&](Subgraph &&subgraph) {
            Cost min_cost = cost(subgraph, m_marked, m_costs);
            subgraphs.emplace_back(min_cost, std::move(subgraph));
            return false;
        });

        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

        if (!subgraphs.empty() && subgraphs[0].first == invalid_cost) {
            state.set_unsolvable();
            return;
        }

        for (auto &[cost, subgraph] : subgraphs) {

            // In the for loop a subgraph can be invalidated when bound_graph is modified.
            // Therefore it has to checked whether subgraph is still valid and does not share a vertex pair with bound_graph.
            const auto vertexPairs = subgraph.vertexPairs();
            bool touches_bound = std::any_of(vertexPairs.begin(), vertexPairs.end(), [&](VertexPair xy) { return !m_marked[xy] && m_bound_graph.has_edge(xy); });

            if (!touches_bound) {
                // std::cout << "after_mark_and_edit insert " << subgraph;
                insert_into_graph(subgraph, m_marked, m_bound_graph);
                state.insert({cost, std::move(subgraph)});
                // std::cout << "\t => " << state << "\n";
            }
        }
    }

private:
    void optimize_bound(State &state, Cost k) {
        std::mt19937 gen(static_cast<unsigned long>(state.cost()) + state.bound.size());
        
        std::shuffle(state.bound.begin(), state.bound.end(), gen);

        bool improvement_found;
        bool bound_changed;
        size_t rounds_no_improvement = 0;
        do {
            improvement_found = false;
            bound_changed = false;

            for (size_t sg_i = 0; sg_i < state.bound.size(); ++sg_i) {
                find_2_improvement(state, sg_i, gen, improvement_found, bound_changed);
                if (state.cost() > k) {
                    return;
                }
            }

            rounds_no_improvement = improvement_found ? 0 : rounds_no_improvement + 1;
        } while (improvement_found || (rounds_no_improvement < 5 && bound_changed));
    }

    void find_one_improvements(State &state, size_t index) {
        bool found_improvement = false;

        const auto &[subgraph_cost, subgraph] = state.bound[index];
        assert(subgraph_cost == cost(subgraph, m_marked, m_costs));

        remove_from_graph(subgraph, m_marked, m_bound_graph);


        Cost max_cost = subgraph_cost;
        Subgraph max_subgraph(subgraph);


        const auto pairs = get_pairs(subgraph, m_marked);

#ifndef NDEBUG
        {
            bool touches = false;
            for (VertexPair uv : subgraph.vertexPairs())
                if (m_bound_graph.has_edge(uv)) touches = true;
            assert(!touches);
        }
#endif

        for (VertexPair uv : pairs) {
            assert(!m_bound_graph.has_edge(uv));

            finder->find_near(uv, m_bound_graph, [&](Subgraph &&neighbor) {
#ifndef NDEBUG
                {
                    auto vp = neighbor.vertexPairs();
                    assert(std::none_of(vp.begin(), vp.end(), [&](VertexPair xy) { return m_bound_graph.has_edge(xy); }));
                }
#endif
                Cost n_cost = cost(neighbor, m_marked, m_costs);
                if (n_cost > max_cost) {
                    found_improvement = true;
                    max_cost = n_cost;
                    max_subgraph = std::move(neighbor);
                }

                return false;
            });

            // prevent subgraphs including uv to be counted twice
            m_bound_graph.set_edge(uv);
        }

        for (VertexPair uv : pairs)
            m_bound_graph.clear_edge(uv);


        insert_into_graph(max_subgraph, m_marked, m_bound_graph);
        state.replace(index, {max_cost, std::move(max_subgraph)});

    }

    /**
     * Find a 2 improvement for subgraph at index.
     *
     * @param state
     * @param index
     * @param gen
     * @param improvement_found
     * @param bound_changed
     * @return
     */
    void find_2_improvement(State &state, size_t index, std::mt19937 &gen, bool &improvement_found,
                            bool &bound_changed) {
        constexpr Cost invalid_max_cost = std::numeric_limits<Cost>::min();
        constexpr size_t invalid_index = std::numeric_limits<size_t>::max();

        const auto &[subgraph_cost, subgraph] = state.bound[index];
        assert(subgraph_cost == cost(subgraph, m_marked, m_costs));
        assert(subgraph_cost != invalid_cost);

        remove_from_graph(subgraph, m_marked, m_bound_graph);

        // candidates are subgraphs which are only adjacent to subgraph but no other subgraph in the lower bound.
        const auto pairs = get_pairs(subgraph, m_marked);

#ifndef NDEBUG
        {
            bool touches = false;
            for (VertexPair uv : subgraph.vertexPairs())
                if (!m_marked[uv] && m_bound_graph.has_edge(uv)) touches = true;
            assert(!touches);
        }
#endif
        auto[candidates, border] = get_candidates(*finder, pairs, m_bound_graph);


        std::vector<Cost> candidate_costs(candidates.size());
        for (size_t i = 0; i < candidates.size(); ++i) {
            candidate_costs[i] = cost(candidates[i], m_marked, m_costs);
        }

        Cost max_subgraphs_cost = subgraph_cost;
        std::vector<size_t> max_subgraphs;

        std::vector<size_t> plateau_candidates;

        // for each candidate check if
        //   1. the candidate has a larger cost than the current maximum cost or
        //   2. a pair of two candidates can both be inserted and their cost is larger than the current maximum cost.
        for (size_t pair_i = 0; pair_i < pairs.size(); ++pair_i) {
            for (size_t a_i = border[pair_i]; a_i < border[pair_i + 1]; ++a_i) {
                const Subgraph &a = candidates[a_i];

                insert_into_graph(a, m_marked, m_bound_graph);

                Cost max_cost_b = invalid_max_cost;
                size_t max_b_i = invalid_index;

                // for each candidate pair (a, b)
                for (size_t pair_j = pair_i + 1; pair_j < pairs.size(); ++pair_j) {
                    for (size_t b_i = border[pair_j]; b_i < border[pair_j + 1]; ++b_i) {
                        const Subgraph &b = candidates[b_i];

                        bool inserted = try_insert_into_graph(b, m_marked, m_bound_graph);

                        if (inserted) {
                            Cost cost_b = candidate_costs[b_i];
                            if (cost_b > max_cost_b) {
                                max_cost_b = cost_b;
                                max_b_i = b_i;
                            }
                            remove_from_graph(b, m_marked, m_bound_graph);
                        }
                    }
                }

                Cost cost_a = candidate_costs[a_i];
                if (max_b_i == invalid_index) {
                    // no partner found
                    if (cost_a > max_subgraphs_cost) {
                        max_subgraphs_cost = cost_a;
                        max_subgraphs = {a_i};
                    } else if (cost_a == max_subgraphs_cost) {
                        // TODO: plateau search
                        plateau_candidates.push_back(a_i);
                    }
                } else {
                    // partner found
                    if (cost_a + max_cost_b > max_subgraphs_cost) {
                        max_subgraphs_cost = cost_a + max_cost_b;
                        max_subgraphs = {a_i, max_b_i};
                    }
                }
                remove_from_graph(a, m_marked, m_bound_graph);
            }
        }

        if (!max_subgraphs.empty()) {
            size_t a_i = max_subgraphs[0];
            improvement_found = true;
            bound_changed = true;

            // better candidates found
            insert_into_graph(candidates[a_i], m_marked, m_bound_graph);
            state.replace(index, {candidate_costs[a_i], std::move(candidates[a_i])});

            if (max_subgraphs.size() == 2) {
                size_t b_i = max_subgraphs[1];

                insert_into_graph(candidates[b_i], m_marked, m_bound_graph);
                state.insert({candidate_costs[b_i], std::move(candidates[b_i])});
            }

        } else {
            // subgraph is the best
            if (!plateau_candidates.empty()) {
                std::shuffle(plateau_candidates.begin(), plateau_candidates.end(), gen);
                auto a_i = plateau_candidates[0];

                insert_into_graph(candidates[a_i], m_marked, m_bound_graph);
                state.replace(index, {candidate_costs[a_i], std::move(candidates[a_i])});
            } else {
                insert_into_graph(subgraph, m_marked, m_bound_graph);
            }
        }
    }

    /**
     * Insert unmarked vertex pairs of "subgraph" into "graph".
     *
     * @param subgraph
     * @param forbidden
     * @param graph
     */
    static void insert_into_graph(const Subgraph &subgraph, const VertexPairMap<bool> &marked, Graph &graph) {
        for (VertexPair uv : subgraph.vertexPairs()) {
            if (!marked[uv]) {
                assert(!graph.has_edge(uv));
                graph.set_edge(uv);
            }
        }
    }

    /**
     * Remove unmarked vertex pairs of "subgraph" from "graph".
     *
     * @param subgraph
     * @param forbidden
     * @param graph
     */
    static void remove_from_graph(const Subgraph &subgraph, const VertexPairMap<bool> &marked, Graph &graph) {
        for (VertexPair uv : subgraph.vertexPairs()) {
            if (!marked[uv]) {
                assert(graph.has_edge(uv));
                graph.clear_edge(uv);
            }
        }
    }

    /**
     * If possible insert all non forbidden vertex pairs of "subgraph" into "graph". If one vertex pair is already in
     * "graph", all previous insertions are inverted.
     * Returns true if all vertex pairs were successfully inserted, false otherwise.
     *
     * @param subgraph
     * @param forbidden
     * @param graph
     * @return
     */
    static bool try_insert_into_graph(const Subgraph &subgraph, const VertexPairMap<bool> &forbidden, Graph &graph) {
        const auto vertexPairs = subgraph.vertexPairs();

        bool touches = std::any_of(vertexPairs.begin(), vertexPairs.end(), [&](VertexPair uv) {
            return !forbidden[uv] && graph.has_edge(uv);
        });


        if (!touches) {
            for (VertexPair uv : subgraph.vertexPairs())
                if (!forbidden[uv])
                    graph.set_edge(uv);

            return true;
        } else {
            return false;
        }
    }

    /**
     * Get candidate forbidden subgraphs adjacent to a pair in "pairs".
     * Border is a separation of the candidates by the pairs.
     *
     * @param finder
     * @param pairs
     * @param bound_graph
     * @return candidates and a border array partitioning the candidates by vertex pairs.
     */
    static std::pair<std::vector<Subgraph>, std::vector<size_t>>
    get_candidates(FinderI &finder, const std::vector<VertexPair> &pairs, Graph &bound_graph) {

        // Precondition: subgraph is removed from bound_graph.
        assert(std::none_of(pairs.begin(), pairs.end(), [&](VertexPair uv) { return bound_graph.has_edge(uv); }));

        std::vector<Subgraph> candidates;
        std::vector<size_t> border(pairs.size() + 1);

        for (size_t i = 0; i < pairs.size(); ++i) {
            VertexPair uv = pairs[i];
            assert(!bound_graph.has_edge(uv));

            finder.find_near(uv, bound_graph, [&](Subgraph &&neighbor) {
#ifndef NDEBUG
                auto vp = neighbor.vertexPairs();
                assert(std::none_of(vp.begin(), vp.end(), [&](VertexPair xy) { return bound_graph.has_edge(xy); }));
#endif
                candidates.push_back(std::move(neighbor));
                return false;
            });
            border[i + 1] = candidates.size();

            // prevent subgraphs including uv to be counted twice
            bound_graph.set_edge(uv);
        }

        // reset bound_graph
        for (VertexPair uv : pairs) {
            assert(bound_graph.has_edge(uv));
            bound_graph.clear_edge(uv);
        }

        return {std::move(candidates), std::move(border)};
    }

    /**
     * Return editable vertex pairs of "subgraph".
     *
     * @param subgraph
     * @param forbidden
     * @return
     */
    static std::vector<VertexPair> get_pairs(const Subgraph &subgraph, const VertexPairMap<bool> &forbidden) {
        std::vector<VertexPair> pairs;
        for (VertexPair uv : subgraph.vertexPairs()) {
            if (!forbidden[uv])
                pairs.push_back(uv);
        }
        return pairs;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_LOCALSEARCHLOWERBOUND_H
