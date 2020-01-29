//
// Created by jonas on 24.09.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_OPTIONS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_OPTIONS_H


#include <iostream>
#include <yaml-cpp/yaml.h>


namespace Options {
    enum class Selector {
        FirstFound, LeastWeight, MostMarkedPairs, MostAdjacentSubgraphs
    };

    std::istream &operator>>(std::istream &in, Selector &selector);

    std::ostream &operator<<(std::ostream &os, Selector selector);

    YAML::Emitter &operator<<(YAML::Emitter &out, Selector selector);


    enum class FSG {
        P3, P4, C4P4, P5, C5P5, P6, C6P6, C4_C5_2K2, C4_C5_P5_Bowtie_Necktie
    };

    std::istream &operator>>(std::istream &in, FSG &fsg);

    std::ostream &operator<<(std::ostream &os, FSG fsg);

    YAML::Emitter &operator<<(YAML::Emitter &out, FSG fsg);


    enum class LB {
        Trivial, Greedy, SortedGreedy, LocalSearch, LPRelaxation, ILSMWISSolver
    };

    std::istream &operator>>(std::istream &in, LB &lower_bound);

    std::ostream &operator<<(std::ostream &os, LB lower_bound);

    YAML::Emitter &operator<<(YAML::Emitter &out, LB lower_bound);


    enum class SolverType {
        FPT, ILP
    };

    std::istream &operator>>(std::istream &in, SolverType &type);

    std::ostream &operator<<(std::ostream &os, SolverType type);

    YAML::Emitter &operator<<(YAML::Emitter &out, SolverType type);


    enum class FPTSearchStrategy {
        Fixed, PrunedDelta, Exponential, IncrementByMinCost, IncrementByMultiplier
    };

    std::istream &operator>>(std::istream &in, FPTSearchStrategy &strategy);

    std::ostream &operator<<(std::ostream &os, FPTSearchStrategy strategy);

    YAML::Emitter &operator<<(YAML::Emitter &out, FPTSearchStrategy strategy);
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_OPTIONS_H
