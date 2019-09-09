//
// Created by jonas on 19.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FPTSOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FPTSOLVER_H


#include "Solver.h"
#include "Configuration.h"
#include "Editor.h"

class FPTSolver : public Solver {

    Cost m_k;
    Options::FSG m_fsg;
public:
    FPTSolver(Cost k, Options::FSG fsg) : m_k(k), m_fsg(fsg) {}

    /**
     * Search for a minimum cost which yields solutions.
     *
     * @param instance
     * @return
     */
    std::vector<Solution> search(Instance instance) {
        auto selector = Options::Selector::FirstEditable;
        auto lower_bound = Options::LB::Greedy;

        Cost k = 0;
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

            Editor editor(instance, selector, m_fsg, lower_bound);

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
    std::vector<Solution> solve(Instance instance) override {
        auto selector = Options::Selector::FirstEditable;
        auto lower_bound = Options::LB::Greedy;

        std::vector<Solution> solutions;


        Editor editor(instance, selector, m_fsg, lower_bound);

        editor.edit(m_k, [&](const std::vector<VertexPair> &edits) {
            solutions.emplace_back(instance, edits);
        }, [](Cost, Cost) {});


        std::sort(solutions.begin(), solutions.end());
        return solutions;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FPTSOLVER_H
