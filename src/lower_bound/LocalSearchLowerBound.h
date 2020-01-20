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

            bool operator==(const Element &other) const { return cost == other.cost && subgraph == other.subgraph; }

            bool operator!=(const Element &other) const { return !(*this == other); }
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
            assert(0 <= m_cost && m_cost < invalid_cost);
            m_cost -= m_bound[index].cost;
            m_bound[index] = std::move(m_bound.back());
            m_bound.pop_back();
            assert(0 <= m_cost && m_cost < invalid_cost);
        }

        /**
         * Replaces the element at the index position with the given element.
         *
         * @param index
         * @param element
         */
        inline void replace(size_t index, Element &&element) {
            assert(0 <= m_cost && m_cost < invalid_cost);
#ifndef NDEBUG
            if (m_bound[index] != element)
                for (const auto &e : m_bound)
                    assert(e != element);
#endif
            m_cost -= m_bound[index].cost;
            m_cost += element.cost;
            m_bound[index] = std::move(element);
            assert(0 <= m_cost && m_cost < invalid_cost);
        }

        inline void insert(Element &&element) {
            assert(0 <= m_cost && m_cost < invalid_cost);
#ifndef NDEBUG
            for (const auto &e : m_bound)
                assert(e != element);
#endif
            m_cost += element.cost;
            m_bound.emplace_back(std::move(element));
            assert(0 <= m_cost && m_cost < invalid_cost);
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
                e.cost = get_subgraph_cost(e.subgraph, marked, costs);
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
    std::mt19937_64 m_gen;
    float m_alpha = 0.7; // percent of plateau candidates being chosen from the ones covering the least other subgraphs instead of all candidates
    size_t m_max_rounds_no_improvement = 5;

    const int verbosity = 0;

public:
    explicit LocalSearchLowerBound(const Instance &instance, const VertexPairMap<bool> &marked,
                                   const SubgraphStats &subgraph_stats, std::shared_ptr<FinderI> finder_ref,
                                   int seed = 0) :
            LowerBoundI(std::move(finder_ref)), m_costs(instance.costs), m_marked(marked),
            m_subgraph_stats(subgraph_stats), m_bound_graph(instance.graph.size()),
            m_gen(static_cast<unsigned long>(seed)) {}

    Cost calculate_lower_bound(Cost k) override;

    Cost get_lower_bound() override { return current_state().cost(); }

    void initialize(Cost /*k*/) override;

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

        // TODO: Check
        initialize_bound_graph(current_state(), m_marked, m_bound_graph);
    }

    State &current_state() {
        assert(!m_states.empty());
        return *m_states.back();
    }

    State &parent_state() {
        assert(m_states.size() > 1);
        return *m_states[m_states.size() - 2];
    }

    void before_mark(VertexPair uv) override;

    void before_edit(VertexPair uv) override;

    void after_edit(VertexPair uv) override;

private:
    static bool bound_graph_is_valid(State &state, const VertexPairMap<bool> &marked, const Graph &bound_graph);

    static bool state_is_valid(State &state, const VertexPairMap<bool> &marked, const VertexPairMap<Cost> &costs);

    static bool bound_is_maximal(FinderI &finder, const Graph &bound_graph);

    static void initialize_bound_graph(const State &state, const VertexPairMap<bool> &marked, Graph &bound_graph);

    static void remove_near_subgraphs_from_bound(State &state, VertexPair uv);

    static void update_near_subgraphs(State &state, VertexPair uv, FinderI &finder,
                                      const VertexPairMap<bool> &marked, const VertexPairMap<Cost> &costs,
                                      Graph &bound_graph);

    static void insert_near_subgraphs_into_bound(State &state, VertexPair uv, FinderI &finder,
                                                 const VertexPairMap<bool> &marked, const VertexPairMap<Cost> &costs,
                                                 Graph &bounded_graph);

    static void insert_subgraphs_into_bound(std::vector<std::pair<Cost, Subgraph>> &&subgraphs,
                                            const VertexPairMap<bool> &marked, State &state, Graph &bound_graph);

    void local_search(State &state, Cost k);

    bool find_one_improvements(State &state, size_t index);

    bool find_two_improvement(State &state, size_t index, bool &bound_changed);

    bool find_omega_improvement(State &state, Cost k);

    static void insert_into_graph(const Subgraph &subgraph, const VertexPairMap<bool> &marked, Graph &graph);

    static void remove_from_graph(const Subgraph &subgraph, const VertexPairMap<bool> &marked, Graph &graph);

    static bool try_insert_into_graph(const Subgraph &subgraph, const VertexPairMap<bool> &marked, Graph &graph);

    static std::pair<std::vector<Subgraph>, std::vector<size_t>> get_candidates(FinderI &finder,
                                                                                const std::vector<VertexPair> &pairs,
                                                                                Graph &bound_graph,
                                                                                const SubgraphStats &subgraph_stats);

    static std::vector<VertexPair> get_pairs(const Subgraph &subgraph, const VertexPairMap<bool> &marked);

    static std::tuple<size_t, size_t> count_neighbors(const SubgraphStats &subgraph_stats,
                                                      const VertexPairMap<bool> &marked, const Subgraph &subgraph);
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_LOCALSEARCHLOWERBOUND_H
