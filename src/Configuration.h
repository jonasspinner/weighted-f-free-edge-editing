//
// Created by jonas on 17.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_CONFIGURATION_H
#define WEIGHTED_F_FREE_EDGE_EDITING_CONFIGURATION_H


#include <string>
#include <iostream>
#include <yaml-cpp/emitter.h>
#include <boost/program_options.hpp>

#include "definitions.h"
#include "options.h"


class Configuration {
public:
    std::string input_path;
    std::string output_path;

    int verbosity = 0;

    // seed and multiplier for input instances
    int seed = 0;
    double multiplier;
    int permutation = 0;

    Options::SolverType solver_type;

    Cost k_max;
    Options::FSG forbidden_subgraphs;

    // FPT parameters
    Options::Selector selector;
    Options::LB lower_bound;
    bool find_all_solutions = true;


    // ILP parameters
    bool sparse_constraints = false;
    bool single_constraints = false;
    int num_threads = 4;
    int timelimit = -1;


    Configuration(Options::FSG forbidden_subgraphs_, std::string input_path_, double multiplier_,
                  Options::SolverType solver_type_, Cost k_max_, Options::Selector selector_, Options::LB lower_bound_) :
            input_path(std::move(input_path_)),
            multiplier(multiplier_),
            solver_type(solver_type_),
            k_max(k_max_),
            forbidden_subgraphs(forbidden_subgraphs_),
            selector(selector_),
            lower_bound(lower_bound_) {}

    boost::program_options::options_description
    options(const std::set<Options::SolverType> &types = {Options::SolverType::FPT, Options::SolverType::ILP});

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const Configuration &config);

    friend std::ostream &operator<<(std::ostream &os, Configuration &config);
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_CONFIGURATION_H
