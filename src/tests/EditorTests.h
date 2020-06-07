//
// Created by jonas on 12.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_EDITORTESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_EDITORTESTS_H


#include <random>

#include "../Configuration.h"


class EditorTests {

    std::mt19937 gen;
    std::string instance_path = "../data/bio/bio-nr-3-size-16.graph";
public:
    explicit EditorTests(int seed = 0) : gen(static_cast<unsigned long>(seed)) {}

    void configurations_have_same_output(Options::FSG fsg, const std::vector<Options::Selector> &selectors,
                                         const std::vector<Options::LB> &lower_bounds, const std::vector<int> &seeds,
                                         double multiplier);

    static void same_output_for_small_zero_cost_instance(const std::vector<Options::Selector> &selectors,
                                                         const std::vector<Options::LB> &lower_bounds,
                                                         const std::vector<int> &seeds);

    void output_is_independent_of_seed(const std::vector<int> &seeds, Options::Selector selector, Options::LB lb);

    void run();
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_EDITORTESTS_H
