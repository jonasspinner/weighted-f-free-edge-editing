//
// Created by jonas on 17.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_CONFIGURATION_H
#define WEIGHTED_F_FREE_EDGE_EDITING_CONFIGURATION_H


#include <string>
#include <iostream>
#include <yaml-cpp/emitter.h>


namespace Options {
    enum class Selector {
        FirstEditable, LeastWeight
    };

    std::istream& operator>>(std::istream& in, Selector& selector);

    std::ostream &operator<<(std::ostream &os, Selector selector);

    YAML::Emitter &operator<<(YAML::Emitter &out, Selector selector);
}


namespace Options {
    enum class FSG {
        P3, P4, P4C4, P5, P5C5, P6, P6C6, C4_C5_2K2, C4_C5_P5_Bowtie_Necktie
    };

    std::istream& operator>>(std::istream& in, FSG& fsg);

    std::ostream &operator<<(std::ostream &os, FSG fsg);

    YAML::Emitter &operator<<(YAML::Emitter &out, FSG fsg);
}


namespace Options {
    enum class LB {
        No, Greedy, LocalSearch, LinearProgram
    };

    std::istream& operator>>(std::istream& in, LB& lower_bound);

    std::ostream &operator<<(std::ostream &os, LB lower_bound);

    YAML::Emitter &operator<<(YAML::Emitter &out, LB lower_bound);
}


class Configuration {
public:
    int k_max;
    Options::Selector selector;
    Options::FSG forbidden;
    Options::LB lower_bound;
    std::string input_path;
    std::string output_path;

    Configuration(int k_max_, Options::Selector selector_, Options::FSG forbidden_, Options::LB lower_bound_,
                  std::string input_path_, std::string output_path_) :
            k_max(k_max_), selector(selector_), forbidden(forbidden_), lower_bound(lower_bound_),
            input_path(std::move(input_path_)), output_path(std::move(output_path_)) {}

    void read_input(int argc, char *argv[]);
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_CONFIGURATION_H
