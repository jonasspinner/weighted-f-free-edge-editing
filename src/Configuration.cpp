#include "Configuration.h"
#include "options.h"


boost::program_options::options_description
Configuration::options(const std::set<Options::SolverType> &types) {
    namespace po = boost::program_options;

    po::options_description cmdline_options;

    po::options_description general("General options");
    general.add_options()
            ("help", "produce help message")
            ("output", po::value<std::string>(&output_path)->default_value(output_path), "output file path")
            ("seed", po::value<int>(&seed)->default_value(seed), "seed for instance permutation")
            ("verbosity", po::value<int>(&verbosity)->default_value(verbosity))
            ("timelimit", po::value<int>(&timelimit)->default_value(timelimit));

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
                ("lower-bound", po::value<Options::LB>(&lower_bound)->default_value(lower_bound))
                ("all", po::value<bool>(&find_all_solutions)->default_value(find_all_solutions))
                ("pre-mark", po::value<bool>(&pre_mark_vertex_pairs)->default_value(pre_mark_vertex_pairs))
                ("search-strategy",
                 po::value<Options::FPTSearchStrategy>(&search_strategy)->default_value(search_strategy));

        cmdline_options.add(fpt_options);
    }

    if (types.count(Options::SolverType::ILP)) {
        po::options_description ilp_options("ILP algorithm options");
        ilp_options.add_options()
                ("constraints", po::value<Options::ILPConstraintGeneration>(&constraints)->default_value(constraints), "if and how constraint generation is restricted")
                ("num-threads", po::value<int>(&num_threads)->default_value(num_threads));

        cmdline_options.add(ilp_options);
    }

    return cmdline_options;
}

YAML::Emitter &operator<<(YAML::Emitter &out, const Configuration &config) {
    using namespace YAML;
    out << BeginMap;
    out << Key << "output" << Value << config.output_path;
    out << Key << "seed" << Value << config.seed;
    out << Key << "solver" << Value << config.solver_type;
    out << Key << "timelimit" << Value << config.timelimit;

    out << Key << "input" << Value << config.input_path;
    out << Key << "permutation" << Value << config.permutation;
    out << Key << "multiplier" << Value << config.multiplier;
    out << Key << "forbidden_subgraphs" << Value << config.forbidden_subgraphs;

    if (config.solver_type == Options::SolverType::FPT) {
        out << Key << "k" << Value << config.k_max;
        out << Key << "selector" << Value << config.selector;
        out << Key << "lower_bound" << Value << config.lower_bound;
        out << Key << "find_all_solutions" << Value << config.find_all_solutions;
        out << Key << "pre_mark_vertex_pairs" << Value << config.pre_mark_vertex_pairs;
        out << Key << "search_strategy" << Value << config.search_strategy;
    } else if (config.solver_type == Options::SolverType::ILP) {
        out << Key << "constraints" << Value << config.constraints;
        out << Key << "num_threads" << Value << config.num_threads;
    }
    out << EndMap;
    return out;
}

std::ostream &operator<<(std::ostream &os, Configuration &config) {
    YAML::Emitter out;
    out << config;
    return os << out.c_str();
}