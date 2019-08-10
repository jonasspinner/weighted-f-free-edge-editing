//
// Created by jonas on 17.07.19.
//

#ifndef CONCEPT_CONFIGURATION_H
#define CONCEPT_CONFIGURATION_H


#include <string>

namespace Options {
    enum Selector {
        FirstEditable, LeastWeight
    };
    enum FSG {
        P3, P4C4
    };
    enum LB {
        No, LocalSearch
    };
}

class Configuration {
public:
    int k_max{100};
    enum SelectorOption {
        FirstEditable, LeastWeight
    } selector;
    enum ForbiddenSubgraphs {
        P3, P4C4
    } forbidden;
    enum LowerBound {
        No, LocalSearch, Greedy
    } lower_bound;
    std::string input_path;
    std::string output_path;

    Configuration(int k_max_, SelectorOption selector_, ForbiddenSubgraphs forbidden_, LowerBound lower_bound_,
                  std::string input_path_, std::string output_path_) :
            k_max(k_max_), selector(selector_), forbidden(forbidden_), lower_bound(lower_bound_),
            input_path(std::move(input_path_)), output_path(std::move(output_path_)) {}

    Configuration(int argc, char *argv[]);
};


#endif //CONCEPT_CONFIGURATION_H
