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

    std::vector<Solution> solve(Instance instance) override {
        auto selector = Options::Selector::FirstEditable;
        auto lower_bound = Options::LB::Greedy;

        std::vector<Solution> solutions;


        Editor editor(instance, selector, m_fsg, lower_bound);

        editor.edit(m_k, [&](const std::vector<VertexPair> &edits) {
            solutions.emplace_back(instance, edits);
        }, [](Cost, Cost) {});


        return solutions;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FPTSOLVER_H
