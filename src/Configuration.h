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
    std::vector<std::string> input_paths;
    std::string output_path;

    int verbosity;

    // seed and multiplier for input instances
    int seed;
    double multiplier;

    Options::SolverType solver_type;

    Cost k_max;
    Options::FSG forbidden_subgraphs;

    // FPT parameters
    Options::Selector selector;
    Options::LB lower_bound;


    Configuration(Cost k_max_, Options::Selector selector_, Options::FSG forbidden_subgraphs_, Options::LB lower_bound_,
                  std::vector<std::string> input_paths_, std::string output_path_, double multiplier_,
                  Options::SolverType solver_type_, int seed_, int verbosity_) :
            input_paths(std::move(input_paths_)),
            output_path(std::move(output_path_)),
            verbosity(verbosity_),
            seed(seed_),
            multiplier(multiplier_),
            solver_type(solver_type_),
            k_max(k_max_),
            forbidden_subgraphs(forbidden_subgraphs_),
            selector(selector_),
            lower_bound(lower_bound_) {}

    void read_input(int argc, char *argv[]) {
        namespace po = boost::program_options;

        auto desc = options();

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            abort();
        }


        if (vm.count("inputs"))
            input_paths = vm["inputs"].as<std::vector<std::string>>();

    }

    boost::program_options::options_description options() {
        namespace po = boost::program_options;

        po::options_description desc("Allowed options");
        desc.add_options()
                ("help", "produce help message")
                ("inputs", po::value<std::vector<std::string>>()->multitoken(), "paths to input instances")
                ("output", po::value<std::string>(&output_path)->default_value(output_path), "output file path")
                ("seed", po::value<int>(&seed)->default_value(seed), "seed for instance permutation")
                ("multiplier", po::value<double>(&multiplier)->default_value(multiplier),
                 "multiplier for discretization of input instances")
                ("solver", po::value<Options::SolverType>(&solver_type)->default_value(solver_type), "type of solver")
                ("k", po::value<Cost>(&k_max)->default_value(k_max), "maximum editing cost")
                ("F", po::value<Options::FSG>(&forbidden_subgraphs)->default_value(forbidden_subgraphs),
                 "forbidden subgraphs")
                ("selector", po::value<Options::Selector>(&selector)->default_value(selector))
                ("lower_bound", po::value<Options::LB>(&lower_bound)->default_value(lower_bound));

        return desc;
    }

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const Configuration &config) {
        using namespace YAML;
        out << BeginMap;
        out << Key << "inputs" << Value << Flow << BeginSeq;
        for (const auto &path : config.input_paths) out << path;
        out << EndSeq;
        out << Key << "output" << Value << config.output_path;
        out << Key << "multiplier" << Value << config.multiplier;
        out << Key << "solver" << Value << config.solver_type;
        out << Key << "k" << Value << config.k_max;
        out << Key << "forbidden_subgraphs" << Value << config.forbidden_subgraphs;
        out << Key << "selector" << Value << config.selector;
        out << Key << "lower_bound" << Value << config.lower_bound;
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
