//
// Created by jonas on 25.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LOCALSEARCHLOWERBOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LOCALSEARCHLOWERBOUND_H


#include <random>
#include "../interfaces/LowerBoundI.h"
#include "../interfaces/FinderI.h"
#include "../consumer/SubgraphStats.h"


class LocalSearchLowerBound : public LowerBoundI {
private:
    class State {
    public:
        struct Element {
            Cost cost;
            Subgraph subgraph;
        };
    private:
        std::vector<Element> m_bound;
        Cost m_cost = 0;
    public:
        [[nodiscard]] inline const std::vector<Element> &bound() const { return m_bound; }

        [[nodiscard]] inline const Element &bound(size_t index) const { return m_bound[index]; }

        [[nodiscard]] inline Cost cost() const { return m_cost; }

        /**
         * Shuffle the order of the subgraphs in the bound.
         *
         * @param g A random number generator.
         */
        template<class URBG>
        inline void shuffle(URBG &&g) {
            std::shuffle(m_bound.begin(), m_bound.end(), g);
        }

        /**
         * Removes the subgraph at the given index. The position is filled by the last subgraph.
         *
         * @param index
         */
        inline void remove(size_t index) {
            m_cost -= m_bound[index].cost;
            m_bound[index] = std::move(m_bound.back());
            m_bound.pop_back();
        }

        /**
         * Replaces the element at the index position with the given element.
         *
         * @param index
         * @param element
         */
        inline void replace(size_t index, Element &&element) {
            m_cost -= m_bound[index].cost;
            m_cost += element.cost;
            m_bound[index] = std::move(element);
        }

        inline void insert(Element &&element) {
            m_cost += element.cost;
            m_bound.emplace_back(std::move(element));
        }

        inline void set_unsolvable() {
            m_bound.clear();
            m_cost = invalid_cost;
        }

        [[nodiscard]] inline bool solvable() const { return m_cost != invalid_cost; }

        /**
         * Recalculate the cost of the bound from the marked vertex pairs and costs.
         *
         * @param marked
         * @param costs
         */
        void recalculate(const VertexPairMap<bool> &marked, const VertexPairMap<Cost> &costs) {
            m_cost = 0;
            for (auto &e : m_bound) {
                e.cost = ::get_subgraph_cost(e.subgraph, marked, costs);
                m_cost += e.cost;
                if (e.cost == invalid_cost) {
                    set_unsolvable();
                    break;
                }
            }

#ifndef NDEBUG
            if (solvable()) {
                Cost sum = 0;
                for (const auto &e : m_bound) {
                    assert(e.cost == ::get_subgraph_cost(e.subgraph, marked, costs));
                    sum += e.cost;
                }
                assert(m_cost == sum);
            }
#endif
        }

        /**
         * Initializes the bound_graph. Every unmarked vertex pair of a subgraph in the bound becomes an edge in the graph.
         *
         * @param marked
         * @param bound_graph
         */
        void initialize_bound_graph(const VertexPairMap<bool> &marked, Graph &bound_graph) {
            bound_graph.clearEdges();

            for (const auto &[cost, subgraph] : m_bound) {
                for (VertexPair uv : subgraph.vertexPairs()) {
                    if (!marked[uv]) {
                        assert(!bound_graph.hasEdge(uv));
                        bound_graph.setEdge(uv);
                    }
                }
            }

#ifndef NDEBUG
            VertexPairMap<bool> debug(bound_graph.size());
            for (const auto&[cost, subgraph] : m_bound)
                for (VertexPair uv : subgraph.vertexPairs())
                    if (!marked[uv]) {
                        assert(bound_graph.hasEdge(uv));
                        assert(!debug[uv]);
                        debug[uv] = true;
                    }

            for (VertexPair uv : bound_graph.vertexPairs())
                if (marked[uv]) {
                    assert(!bound_graph.hasEdge(uv));
                    assert(!debug[uv]);
                } else {
                    assert(bound_graph.hasEdge(uv) == static_cast<bool>(debug[uv]));
                }
#endif
        }

        void
        assert_valid(const VertexPairMap<bool> &marked, const Graph &bound_graph, const VertexPairMap<Cost> &costs) {
#ifndef NDEBUG
            if (!solvable()) return;
            // every subgraph is valid
            //for (const auto &[_, subgraph] : bound) {
            //    verify(subgraph, graph);
            //}

            // costs are matching
            Cost sum = 0;
            for (const auto &e : m_bound) {
                assert(e.cost == get_subgraph_cost(e.subgraph, marked, costs));
                sum += e.cost;
            }
            assert(m_cost == sum);

            VertexPairMap<bool> debug(bound_graph.size());

            for (const auto&[cost, subgraph] : m_bound) {
                for (VertexPair uv : subgraph.vertexPairs()) {
                    if (!marked[uv]) {
                        assert(bound_graph.hasEdge(uv));
                        assert(!debug[uv]);
                        debug[uv] = true;
                    }
                }
            }

            for (VertexPair uv : bound_graph.vertexPairs()) {
                if (marked[uv]) {
                    assert(!bound_graph.hasEdge(uv));
                    assert(!debug[uv]);
                } else {
                    assert(bound_graph.hasEdge(uv) == static_cast<bool>(debug[uv]));
                }
            }
#endif
        }

        friend std::ostream &operator<<(std::ostream &os, const State &state) {
            if (!state.solvable()) return os << "unsolvable";
            std::cout << state.cost() << ":";
            for (const auto &[cost, subgraph] : state.bound())
                std::cout << " (" << subgraph << ", " << cost << ")";
            return os;
        }
    };

    const VertexPairMap<Cost> &m_costs;
    const VertexPairMap<bool> &m_marked;
    const SubgraphStats &m_subgraph_stats;

    Graph m_bound_graph;
    std::vector<std::unique_ptr<State>> m_states;
    std::mt19937_64 gen;

    const int verbosity = 0;

public:
    explicit LocalSearchLowerBound(const Instance &instance, const VertexPairMap<bool> &marked,
                                   const SubgraphStats &subgraph_stats, std::shared_ptr<FinderI> finder_ref,
                                   int seed = 0) :
            LowerBoundI(std::move(finder_ref)), m_costs(instance.costs), m_marked(marked),
            m_subgraph_stats(subgraph_stats), m_bound_graph(instance.graph.size()), gen(seed) {}

    /**
     * Returns a lower bound on the editing cost.
     *
     * @param k Remaining editing cost.
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
            local_search(state, k);
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
    void initialize(Cost /*k*/) override {
        using Element = State::Element;
        m_states.push_back(std::make_unique<State>());
        State &state = *m_states.back();

        m_bound_graph.clearEdges();


        std::vector<Element> subgraphs;

        finder->find([&](Subgraph &&subgraph) {
            Cost min_cost = get_subgraph_cost(subgraph, m_marked, m_costs);
            subgraphs.push_back({min_cost, std::move(subgraph)});
            return false;
        });

        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const Element &lhs, const Element &rhs) { return lhs.cost > rhs.cost; });

        if (!subgraphs.empty() && subgraphs[0].cost == invalid_cost) {
            state.set_unsolvable();
            return;
        }

        for (auto&[cost, subgraph] : subgraphs) {
            bool inserted = try_insert_into_graph(subgraph, m_marked, m_bound_graph);

            if (inserted) {
                assert(cost != invalid_cost);
                state.insert({cost, std::move(subgraph)});
            }
        }

        /*
        auto unsolvable = finder->find([&](Subgraph &&subgraph) {
            Cost min_cost = cost(subgraph, m_marked, m_costs);
            bool inserted = try_insert_into_graph(subgraph, m_marked, m_bound_graph);

            if (min_cost == invalid_cost)
                return true;

            if (inserted) {
                assert(min_cost != invalid_cost);
                state.insert({min_cost, std::move(subgraph)});
            }

            return false;
        });

        if (unsolvable) {
            state.set_unsolvable();
            return;
        }*/

        if (verbosity) std::cout << "[" << state.cost() << "] greedy lb\n";
        local_search(state, std::numeric_limits<Cost>::max());
        if (verbosity) std::cout << "[" << state.cost() << "] after first local search\n";
    }

    /**
     * Copy the state from the previous level.
     *
     * @param k
     * @return
     */
    void push_state(Cost /*k*/) override {
        assert(!m_states.empty());
        m_states.push_back(std::make_unique<State>(*m_states.back()));
    }

    void pop_state() override {
        assert(!m_states.empty());
        m_states.pop_back();
    }

    State &current_state() {
        assert(!m_states.empty());
        return *m_states.back();
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

        for (size_t i = 0; i < state.bound().size();) {
            const auto &[cost, subgraph] = state.bound(i);
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

    /**
     * Find forbidden subgraphs near uv and insert them. Subgraphs with higher minimum editing cost are preferred.
     *
     * @param uv
     */
    void after_mark_and_edit(VertexPair uv) override {
        auto &state = current_state();
        assert(state.solvable());

        // TODO: Check if needed
        state.initialize_bound_graph(m_marked, m_bound_graph);

        std::vector<std::pair<Cost, Subgraph>> subgraphs;

        // The finder iterates over subgraphs having u and v as vertices.
        finder->find_near(uv, m_bound_graph, [&](Subgraph &&subgraph) {
            Cost min_cost = get_subgraph_cost(subgraph, m_marked, m_costs);
            subgraphs.emplace_back(min_cost, std::move(subgraph));
            return false;
        });

        std::sort(subgraphs.begin(), subgraphs.end(),
                  [](const auto &lhs, const auto &rhs) { return lhs.first > rhs.first; });

        // if one forbidden subgraph has no unmarked pairs, the instance is unsolvable.
        if (!subgraphs.empty() && subgraphs[0].first == invalid_cost) {
            state.set_unsolvable();
            return;
        }

        for (auto &[cost, subgraph] : subgraphs) {

            // In the for loop a subgraph can be invalidated when bound_graph is modified.
            // Therefore it has to checked whether subgraph is still valid and does not share a vertex pair with bound_graph.
            const auto vertexPairs = subgraph.vertexPairs();
            bool touches_bound = std::any_of(vertexPairs.begin(), vertexPairs.end(), [&](VertexPair xy) {
                return !m_marked[xy] && m_bound_graph.hasEdge(xy);
            });

            if (!touches_bound) {
                // std::cout << "after_mark_and_edit insert " << subgraph;
                insert_into_graph(subgraph, m_marked, m_bound_graph);
                state.insert({cost, std::move(subgraph)});
                // std::cout << "\t => " << state << "\n";
            }
        }
    }

private:
    /**
     * Tries to optimise the bound. Stops if the lower bound is larger than k or the bound does not improve in five consecutive rounds.
     *
     * @param state
     * @param k
     */
    void local_search(State &state, Cost k) {
        state.shuffle(gen);

        bool improvement_found, bound_changed;
        size_t rounds_no_improvement = 0;

        do {
            // a single round consists of a loop over all subgraphs in the bound and trying to find an improvement depending on the current mode.
            improvement_found = bound_changed = false;

            if (verbosity) std::cout << "[" << state.cost() << "] round mode=0\n";

            for (size_t index = 0; state.cost() <= k && index < state.bound().size(); ++index)
                improvement_found |= find_one_improvements(state, index);

            if (!improvement_found) {
                if (verbosity) std::cout << "[" << state.cost() << "] round mode=1\n";
                for (size_t index = 0; state.cost() <= k && index < state.bound().size(); ++index)
                    improvement_found |= find_two_improvement(state, index, bound_changed);
            }

            if (!improvement_found) {
                if (verbosity) std::cout << "[" << state.cost() << "] round mode=2\n";
                improvement_found = find_omega_improvement(state, k);
            }

            rounds_no_improvement = improvement_found ? 0 : rounds_no_improvement + 1;
            // improvement_found || (rounds_no_improvement < 5 && bound_changed)
        } while (state.cost() <= k && (improvement_found || (rounds_no_improvement < 5 && bound_changed)));
        if (verbosity) std::cout << "[" << state.cost() << "] end local search\n";
    }

    /**
     * Try to find a (1, 1) swap for the subgraph at the index which improves the lower bound.
     *
     * For the subgraph currently in the lower bound at the given index do the following. List other forbidden subgraphs
     * which are adjacent to the subgraph but not to an other subgraph in the lower bound as candidates.
     * + If a candidate has a larger cost than the original subgraph, choose it instead.
     *
     * Complexity: O((p over 2) * find_near + |candidates|)
     *
     * @param state
     * @param index
     */
    bool find_one_improvements(State &state, size_t index) {
        bool found_improvement = false;

        const auto &[subgraph_cost, subgraph] = state.bound(index);
        assert(subgraph_cost == get_subgraph_cost(subgraph, m_marked, m_costs));

        remove_from_graph(subgraph, m_marked, m_bound_graph);


        Cost max_cost = subgraph_cost;
        Subgraph max_subgraph(subgraph);


        const auto pairs = get_pairs(subgraph, m_marked);

#ifndef NDEBUG
        {
            bool touches = false;
            for (VertexPair uv : subgraph.vertexPairs())
                if (m_bound_graph.hasEdge(uv)) touches = true;
            assert(!touches);
        }
#endif

        for (VertexPair uv : pairs) {
            assert(!m_bound_graph.hasEdge(uv));

            finder->find_near(uv, m_bound_graph, [&](Subgraph &&neighbor) {
#ifndef NDEBUG
                {
                    auto vp = neighbor.vertexPairs();
                    assert(std::none_of(vp.begin(), vp.end(),
                                        [&](VertexPair xy) { return m_bound_graph.hasEdge(xy); }));
                }
#endif
                Cost n_cost = get_subgraph_cost(neighbor, m_marked, m_costs);
                if (n_cost > max_cost) {
                    found_improvement = true;
                    max_cost = n_cost;
                    max_subgraph = std::move(neighbor);
                }

                return false;
            });

            // prevent subgraphs including uv to be counted twice
            m_bound_graph.setEdge(uv);
        }

        for (VertexPair uv : pairs)
            m_bound_graph.clearEdge(uv);

        if (max_cost > subgraph_cost && verbosity > 1)
            std::cout << "found (1, 1) swap " << std::setw(4) << max_cost - subgraph_cost << ", " << subgraph_cost
                      << " => " << max_cost << ", " << subgraph << " => " << max_subgraph << "\n";

        insert_into_graph(max_subgraph, m_marked, m_bound_graph);
        state.replace(index, {max_cost, std::move(max_subgraph)});

        return found_improvement;
    }

    /**
     * Try to find a (1, 2) swap for the subgraph at the index which improves the lower bound.
     *
     * For the subgraph currently in the lower bound at the given index do the following. List other forbidden subgraphs
     * which are adjacent to the subgraph but not to an other subgraph in the lower bound as candidates. Test for each
     * pair of candidates whether they are adjacent.
     * + From all pairs of candidates which are not adjacent, choose the pair with the largest combined cost.
     * + If single candidate has a larger cost choose that.
     * + If no improvement can be found choose a candidate with the same cost as the original subgraph at random.
     * + If no candidate can be found, let the original subgraph remain.
     *
     * Complexity: O((p over 2) * find_near + |candidates|^2)
     *
     * @param state
     * @param index
     * @param gen
     * @param bound_changed
     * @return
     */
    bool find_two_improvement(State &state, size_t index, bool &bound_changed) {
        constexpr Cost invalid_max_cost = std::numeric_limits<Cost>::min();
        constexpr size_t invalid_index = std::numeric_limits<size_t>::max();

        bool found_improvement = false;

        const auto &[subgraph_cost, subgraph] = state.bound(index);
        assert(subgraph_cost == get_subgraph_cost(subgraph, m_marked, m_costs));
        assert(subgraph_cost != invalid_cost);

        remove_from_graph(subgraph, m_marked, m_bound_graph);

        // candidates are subgraphs which are only adjacent to subgraph but no other subgraph in the lower bound.
        const auto pairs = get_pairs(subgraph, m_marked);

#ifndef NDEBUG
        {
            bool touches = false;
            for (VertexPair uv : subgraph.vertexPairs())
                if (!m_marked[uv] && m_bound_graph.hasEdge(uv)) touches = true;
            assert(!touches);
        }
#endif
        auto[candidates, border] = get_candidates(*finder, pairs, m_bound_graph);


        std::vector<Cost> candidate_costs(candidates.size());
        for (size_t i = 0; i < candidates.size(); ++i) {
            candidate_costs[i] = get_subgraph_cost(candidates[i], m_marked, m_costs);
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
            found_improvement = true;
            bound_changed = true;

            if (verbosity > 1)
                std::cout << "found (1, " << max_subgraphs.size() << ") swap " << std::setw(4)
                          << max_subgraphs_cost - subgraph_cost << ", " << subgraph_cost << " => " << max_subgraphs_cost
                          << ", " << subgraph << " => ...\n";

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
                bound_changed = true;

                std::shuffle(plateau_candidates.begin(), plateau_candidates.end(), gen);
                auto a_i = plateau_candidates[0];

                if (verbosity > 1)
                    std::cout << "made (1, 1) swap for plateau search " << std::setw(4) << 0 << ", " << subgraph
                              << " => " << candidates[a_i] << "\n";

                insert_into_graph(candidates[a_i], m_marked, m_bound_graph);
                state.replace(index, {candidate_costs[a_i], std::move(candidates[a_i])});
            } else {
                insert_into_graph(subgraph, m_marked, m_bound_graph);
            }
        }
        return found_improvement;
    }

    /**
     * Try to find a (\omega, 1) swap which improves the lower bound.
     *
     * Find all forbidden subgraph in the graph. Check which subgraphs in the lower bound are adjacent to it. If the
     * cost of the single subgraph is larger than the sum of the adjacent subgraphs replace them.
     *
     * Complexity: O(find * (p choose 2))
     *
     * @param state
     * @param k
     * @return
     */
    bool find_omega_improvement(State &state, Cost k) {
        constexpr size_t invalid_index = std::numeric_limits<size_t>::max();

        bool found_improvement = false;

        VertexPairMap<bool> used(m_bound_graph.size());
        VertexPairMap<size_t> subgraph_index(m_bound_graph.size(), invalid_index);

        auto index_assign = [&](const Subgraph &subgraph, size_t i) {
            for (VertexPair uv : subgraph.vertexPairs())
                if (!m_marked[uv])
                    subgraph_index[uv] = i;
        };

        auto candidate_neighbors_indices = [&](const Subgraph &subgraph) {
            std::vector<size_t> is;
            auto vp = subgraph.vertexPairs();
            std::transform(vp.begin(), vp.end(), std::back_inserter(is),
                           [&](VertexPair uv) { return subgraph_index[uv]; });
            is.erase(std::remove_if(is.begin(), is.end(), [](size_t i) { return i == invalid_index; }), is.end());
            std::sort(is.begin(), is.end());
            is.erase(std::unique(is.begin(), is.end()), is.end());
            return is;
        };

        for (size_t i = 0; i < state.bound().size(); ++i)
            index_assign(state.bound(i).subgraph, i);


        finder->find([&](Subgraph &&subgraph) {
            auto subgraph_cost = get_subgraph_cost(subgraph, m_marked, m_costs);
            assert(subgraph_cost != invalid_cost);

            for (VertexPair uv : subgraph.vertexPairs())
                if (!m_marked[uv])
                    used[uv] = true;

            Cost sum = 0;
            size_t count = 0;

            auto is = candidate_neighbors_indices(subgraph);
            for (auto i : is) {
                sum += state.bound(i).cost;
                count++;
            }

            if (subgraph_cost > sum) {
                found_improvement = true;
                if (verbosity > 1)
                    std::cout << "found (" << count << ", 1) swap " << std::setw(4) << subgraph_cost - sum << ", "
                              << sum << " => " << subgraph_cost << ", ... => " << subgraph << "\n";

                for (VertexPair uv : subgraph.vertexPairs()) {
                    auto i = subgraph_index[uv];
                    if (i != invalid_index) {
                        const auto &bound_subgraph = state.bound(i).subgraph;

                        index_assign(bound_subgraph, invalid_index);
                        remove_from_graph(bound_subgraph, m_marked, m_bound_graph);
                        state.remove(i);
                        index_assign(state.bound(i).subgraph, i);
                    }
                }

                for (VertexPair uv : subgraph.vertexPairs())
                    used[uv] = false;

                index_assign(subgraph, state.bound().size());
                insert_into_graph(subgraph, m_marked, m_bound_graph);
                state.insert({subgraph_cost, std::move(subgraph)});
            } else {
                for (VertexPair uv : subgraph.vertexPairs())
                    used[uv] = false;
            }
            return state.cost() > k;
        });

        return found_improvement;
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
                assert(!graph.hasEdge(uv));
                graph.setEdge(uv);
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
                assert(graph.hasEdge(uv));
                graph.clearEdge(uv);
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
            return !forbidden[uv] && graph.hasEdge(uv);
        });


        if (!touches) {
            for (VertexPair uv : subgraph.vertexPairs())
                if (!forbidden[uv])
                    graph.setEdge(uv);

            return true;
        }
        return false;
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
        assert(std::none_of(pairs.begin(), pairs.end(), [&](VertexPair uv) { return bound_graph.hasEdge(uv); }));

        std::vector<Subgraph> candidates;
        std::vector<size_t> border(pairs.size() + 1);

        for (size_t i = 0; i < pairs.size(); ++i) {
            VertexPair uv = pairs[i];
            assert(!bound_graph.hasEdge(uv));

            finder.find_near(uv, bound_graph, [&](Subgraph &&neighbor) {
#ifndef NDEBUG
                auto vp = neighbor.vertexPairs();
                assert(std::none_of(vp.begin(), vp.end(), [&](VertexPair xy) { return bound_graph.hasEdge(xy); }));
#endif
                candidates.push_back(std::move(neighbor));
                return false;
            });
            border[i + 1] = candidates.size();

            // prevent subgraphs including uv to be counted twice
            bound_graph.setEdge(uv);
        }

        // reset bound_graph
        for (VertexPair uv : pairs) {
            assert(bound_graph.hasEdge(uv));
            bound_graph.clearEdge(uv);
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
