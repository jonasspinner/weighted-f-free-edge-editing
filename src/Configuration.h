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


namespace Options {
    enum class Selector {
        FirstEditable, LeastWeight, MostMarked
    };

    std::istream &operator>>(std::istream &in, Selector &selector);

    std::ostream &operator<<(std::ostream &os, Selector selector);

    YAML::Emitter &operator<<(YAML::Emitter &out, Selector selector);


    enum class FSG {
        P3, P4, P4C4, P5, P5C5, P6, P6C6, C4_C5_2K2, C4_C5_P5_Bowtie_Necktie
    };

    std::istream &operator>>(std::istream &in, FSG &fsg);

    std::ostream &operator<<(std::ostream &os, FSG fsg);

    YAML::Emitter &operator<<(YAML::Emitter &out, FSG fsg);


    enum class LB {
        No, Greedy, LocalSearch, LinearProgram
    };

    std::istream &operator>>(std::istream &in, LB &lower_bound);

    std::ostream &operator<<(std::ostream &os, LB lower_bound);

    YAML::Emitter &operator<<(YAML::Emitter &out, LB lower_bound);


    enum class SolverType {
        FPT, ILP
    };

    std::istream &operator>>(std::istream &os, SolverType &type);

    std::ostream &operator<<(std::ostream &os, SolverType type);

    YAML::Emitter &operator<<(YAML::Emitter &out, SolverType type);
}


class Configuration {
public:
    std::string input_path;
    std::string output_path;

    int verbosity;

    // seed and multiplier for input instances
    int seed;
    double multiplier;
    int permutation;

    Options::SolverType solver_type;

    Cost k_max;
    Options::FSG forbidden_subgraphs;

    // FPT parameters
    Options::Selector selector;
    Options::LB lower_bound;
    bool find_all_solutions;


    Configuration(Cost k_max_, Options::Selector selector_, Options::FSG forbidden_subgraphs_, Options::LB lower_bound_,
                  std::string input_path_, std::string output_path_, double multiplier_,
                  Options::SolverType solver_type_, int seed_, int verbosity_, bool find_all_solutions_) :
            input_path(std::move(input_path_)),
            output_path(std::move(output_path_)),
            verbosity(verbosity_),
            seed(seed_),
            multiplier(multiplier_),
            permutation(0),
            solver_type(solver_type_),
            k_max(k_max_),
            forbidden_subgraphs(forbidden_subgraphs_),
            selector(selector_),
            lower_bound(lower_bound_),
            find_all_solutions(find_all_solutions_) {}

    boost::program_options::options_description
    options(const std::set<Options::SolverType> &types = {Options::SolverType::FPT, Options::SolverType::ILP}) {
        namespace po = boost::program_options;

        po::options_description cmdline_options;

        po::options_description general("General options");
        general.add_options()
                ("help", "produce help message")
                ("output", po::value<std::string>(&output_path)->default_value(output_path), "output file path")
                ("seed", po::value<int>(&seed)->default_value(seed), "seed for instance permutation");

        if (types.size() > 1) {
            general.add_options()
                    ("solver", po::value<Options::SolverType>(&solver_type)->default_value(solver_type),
                     "type of solver");
        }

        cmdline_options.add(general);


        po::options_description problem_options("Problem options");
        problem_options.add_options()
                ("input", po::value<std::string>(&input_path)->default_value(input_path), "path to input instance")
                ("permutation", po::value<int>(&permutation)->default_value(permutation), "seed for input permutation")
                ("multiplier", po::value<double>(&multiplier)->default_value(multiplier),
                 "multiplier for discretization of input instances")
                ("F", po::value<Options::FSG>(&forbidden_subgraphs)->default_value(forbidden_subgraphs),
                 "forbidden subgraphs");

        cmdline_options.add(problem_options);


        if (types.count(Options::SolverType::FPT)) {
            po::options_description fpt_options("FPT algorithm options");
            fpt_options.add_options()
                    ("k", po::value<Cost>(&k_max)->default_value(k_max), "maximum editing cost")
                    ("selector", po::value<Options::Selector>(&selector)->default_value(selector))
                    ("lower_bound", po::value<Options::LB>(&lower_bound)->default_value(lower_bound));

            cmdline_options.add(fpt_options);
        }

        if (types.count(Options::SolverType::ILP)) {
            po::options_description ilp_options("ILP algorithm options");
            ilp_options.add_options();

            cmdline_options.add(ilp_options);
        }

        return cmdline_options;
    }

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const Configuration &config) {
        using namespace YAML;
        out << BeginMap;
        out << Key << "output" << Value << config.output_path;
        out << Key << "seed" << Value << config.seed;
        out << Key << "solver" << Value << config.solver_type;

        out << Key << "input" << Value << config.input_path;
        out << Key << "permutation" << Value << config.permutation;
        out << Key << "multiplier" << Value << config.multiplier;
        out << Key << "forbidden_subgraphs" << Value << config.forbidden_subgraphs;

        if (config.solver_type == Options::SolverType::FPT) {
            out << Key << "k" << Value << config.k_max;
            out << Key << "selector" << Value << config.selector;
            out << Key << "lower_bound" << Value << config.lower_bound;
        }
        out << EndMap;
        return out;
    }

    friend std::ostream &operator<<(std::ostream &os, Configuration &config) {
        YAML::Emitter out;
        out << config;
        return os << out.c_str();
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_CONFIGURATION_H
