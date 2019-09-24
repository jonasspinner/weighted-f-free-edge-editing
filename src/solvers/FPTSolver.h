//
// Created by jonas on 19.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FPTSOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FPTSOLVER_H


#include "Solver.h"
#include "../Configuration.h"
#include "../Editor.h"

class FPTSolver : public Solver {

    Configuration m_config;
    Statistics m_stats;
public:
    explicit FPTSolver(Configuration config) : m_config(std::move(config)), m_stats() {}

    /**
     * Search for a minimum cost which yields solutions.
     *
     * @param instance
     * @return
     */
    std::vector<Solution> search(Instance instance) {

        Editor editor(instance, m_config);


        Cost k = editor.initial_lower_bound();
        bool solved;
        std::vector<Solution> solutions;

        auto append_solution = [&](const std::vector<VertexPair> &edits) {
            solutions.emplace_back(instance, edits);
        };

        do {
            Cost min_delta = std::numeric_limits<Cost>::max();
            std::vector<Cost> deltas;
            std::cout << "edit(" << std::setw(10) << k << "):";

            auto handle_delta = [&](Cost pruned_k, Cost lb) { min_delta = std::min(min_delta, lb - pruned_k); deltas.push_back(lb - pruned_k); };

            solved = editor.edit(k, append_solution, handle_delta);

            std::sort(deltas.begin(), deltas.end(), [](Cost lhs, Cost rhs) { return lhs < rhs; });
            Cost delta10 = deltas[std::min(deltas.size() - 1, deltas.size() / 5)];

            std::cout << " min_delta = " << std::setw(10) << min_delta << ", deltas size = " << std::setw(10) << deltas.size() << ", delta10 = " << std::setw(10) << delta10 << "\n";

            if (!solved) k += delta10;
        } while (!solved);

        return solutions;
    }

    /**
     * Solves the given instance.
     *
     * @param instance
     * @return
     */
    Result solve(Instance instance) override {
        std::vector<Solution> solutions;

        Editor editor(instance, m_config);

        auto k_min = editor.initial_lower_bound();
        if (k_min > m_config.k_max)
            return Result::Unsolved();

        bool solved = editor.edit(m_config.k_max, [&](const std::vector<VertexPair> &edits) {
            solutions.emplace_back(instance, edits);
        }, [](Cost, Cost) {});

        m_stats = Statistics(editor.stats());

        if (solved) {
            std::sort(solutions.begin(), solutions.end());
            Solution::filter_inclusion_minimal(solutions);
            return Result::Solved(solutions);
        } else {
            return Result::Unsolved();
        }
    }

    const Statistics &stats() {
        return m_stats;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FPTSOLVER_H
