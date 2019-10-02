//
// Created by jonas on 18.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_MOSTADJACENTSUBGRAPHS_ADAPTED_H
#define WEIGHTED_F_FREE_EDGE_EDITING_MOSTADJACENTSUBGRAPHS_ADAPTED_H


#include "../consumer/SubgraphStats.h"

class MostAdjacentSubgraphs_adapted : public SelectorI {
private:
    struct State {
        std::vector<Subgraph> one_left_subgraphs;
        bool impossible_to_solve;

        State() : impossible_to_solve(false) {}
    };

    struct forbidden_count {
        VertexPair node_pair;
        size_t num_forbidden;

        forbidden_count(VertexPair pair_, size_t num_forbidden_) : node_pair(pair_), num_forbidden(num_forbidden_) {}

        bool operator<(const forbidden_count &other) const {
            return this->num_forbidden > other.num_forbidden;
        }

        operator VertexPair() const {
            return node_pair;
        }
    };

    const VertexPairMap<bool> &m_marked;
    const SubgraphStats &m_subgraph_stats;

    std::vector<std::unique_ptr<State>> m_states;

public:
    /**
     * This class is adapted from the fpt-editing repository.
     */
    MostAdjacentSubgraphs_adapted(std::shared_ptr <FinderI> finder_ptr, const VertexPairMap<bool> &marked, const SubgraphStats &subgraph_stats) : SelectorI(
            std::move(finder_ptr)), m_marked(marked), m_subgraph_stats(subgraph_stats) {
        m_states.push_back(std::make_unique<State>());
    }

private:
    State &current_state() {
        return *m_states.back();
    }

    State &parent_state() {
        return *m_states[m_states.size() - 2];
    }
public:

    void push_state(Cost /*k*/) override {
        m_states.push_back(std::make_unique<State>(*m_states.back()));
    }

    void pop_state() override {
        m_states.pop_back();
    }

    void initialize(Cost /*k*/) override {
        m_states.push_back(std::make_unique<State>());
        State &state = *m_states.back();

        finder->find([&](const Subgraph &path) {
            size_t free = 0;

            for (VertexPair uv : path.vertexPairs())
                if (!m_marked[uv])
                    ++free;


            if (free == 1) {
                state.one_left_subgraphs.push_back(path);
            } else if (free == 0) {
                state.impossible_to_solve = true;
                return true;
            }

            return false;
        });
    }

    void after_edit(VertexPair uv) override {
        State &state = current_state();

        finder->find_near(uv, [&](Subgraph &&subgraph) {
            size_t free = 0;

            for (VertexPair xy : subgraph.vertexPairs())
                if (!m_marked[xy])
                    ++free;

            if (free == 1) {
                state.one_left_subgraphs.push_back(std::move(subgraph));
            } else if (free == 0) {
                state.impossible_to_solve = true;
            }

            return false;
        });
    }

    void after_mark(VertexPair uv) override {
        State &state = parent_state();
        assert(m_marked[uv]);

        finder->find_near(uv, [&](Subgraph &&subgraph) {
            size_t free = 0;

            for (VertexPair xy : subgraph.vertexPairs())
                if (!m_marked[xy])
                    ++free;

            if (free == 1) {
                state.one_left_subgraphs.push_back(std::move(subgraph));
            } else if (free == 0) {
                state.impossible_to_solve = true;
            }

            return false;
        });
    }

    Problem select_problem(Cost k) override {
        State &state = current_state();

        Problem problem;
        problem.solved = (m_subgraph_stats.subgraphCount() == 0);

        if (!problem.solved && k > 0 && !state.impossible_to_solve) {
            while (!state.one_left_subgraphs.empty()) {
                Subgraph subgraph = std::move(state.one_left_subgraphs.back());
                state.one_left_subgraphs.pop_back();

                VertexPair least_subgraphs_pair{0, 1};
                size_t min_subgraph_count = std::numeric_limits<size_t>::max();
                for (VertexPair uv : subgraph.vertexPairs()) {
                    auto subgraph_count = m_subgraph_stats.subgraphCount(uv);
                    if (subgraph_count < min_subgraph_count) {
                        min_subgraph_count = subgraph_count;
                        least_subgraphs_pair = uv;
                    }
                }

                bool is_valid_subgraph = finder->find_near(least_subgraphs_pair, [&](Subgraph &&other) {
                    return subgraph == other;
                });

                if (is_valid_subgraph) {

                    for (VertexPair uv : subgraph.vertexPairs())
                        if (!m_marked[uv])
                            problem.pairs.emplace_back(uv);

                    assert(problem.pairs.size() == 1);

                    return problem;
                }
            }

            size_t max_subgraphs = 0;
            std::vector<VertexPair> pairs;

            for (VertexPair uv : Graph::VertexPairs(m_marked.size())) {
                auto subgraph_count = m_subgraph_stats.subgraphCount(uv);
                if (subgraph_count > max_subgraphs) {
                    max_subgraphs = subgraph_count;
                    pairs.clear();
                }
                if (subgraph_count == max_subgraphs && (subgraph_count > 1 || pairs.empty())) {
                    pairs.emplace_back(uv);
                }
            }

            std::vector<forbidden_count> best_pairs, current_pairs;
            for (VertexPair uv : pairs) {
                // TODO: restrict listing to exclude already listed subgraphs?
                finder->find_near(uv, [&](const Subgraph &subgraph) {
                    current_pairs.clear();

                    for (VertexPair xy : subgraph.vertexPairs())
                        if (!m_marked[xy])
                            current_pairs.emplace_back(xy, m_subgraph_stats.subgraphCount(xy));

                    // assert(current_pairs.size() > 1);

                    std::sort(current_pairs.begin(), current_pairs.end());

                    if (best_pairs.empty()) {
                        best_pairs = std::move(current_pairs);
                    } else {
                        size_t bi = 0, ci = 0;

                        while (bi + 1 < best_pairs.size() && ci + 1 < current_pairs.size() &&
                               best_pairs[bi].num_forbidden == current_pairs[ci].num_forbidden) {
                            ++bi;
                            ++ci;
                        }

                        if (ci + 1 == current_pairs.size() || (bi + 1 != best_pairs.size() &&
                                                               best_pairs[bi].num_forbidden <
                                                               current_pairs[ci].num_forbidden)) {
                            best_pairs = std::move(current_pairs);
                        }
                    }

                    return false;
                });
            }

            for (const auto &pair_count : best_pairs) {
                problem.pairs.emplace_back(pair_count.node_pair);
            }
        }

        return problem;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_MOSTADJACENTSUBGRAPHS_ADAPTED_H
