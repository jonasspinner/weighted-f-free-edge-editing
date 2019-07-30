//
// Created by jonas on 17.07.19.
//

#ifndef CONCEPT_CONFIGURATION_H
#define CONCEPT_CONFIGURATION_H


#include <string>

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
        No, IteratedLocalSearch
    } lower_bound;
    std::string input_path;
    std::string output_path;

    Configuration(int k_max, SelectorOption selector, ForbiddenSubgraphs forbidden, LowerBound lower_bound, std::string input_path,
                  std::string output_path) : k_max(k_max), selector(selector), forbidden(forbidden), lower_bound(lower_bound),
                                             input_path(std::move(input_path)), output_path(std::move(output_path)) {}

    Configuration(int argc, char *argv[]);
};


#endif //CONCEPT_CONFIGURATION_H
