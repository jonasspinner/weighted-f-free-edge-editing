//
// Created by jonas on 17.07.19.
//

#ifndef CONCEPT_CONFIGURATION_H
#define CONCEPT_CONFIGURATION_H


#include <string>

namespace Options {
    enum class Selector {
        FirstEditable, LeastWeight
    };
    enum class FSG {
        P3, P4C4
    };
    enum class LB {
        No, Greedy, LocalSearch
    };
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


#endif //CONCEPT_CONFIGURATION_H
