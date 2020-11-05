#include <iostream>
#include <chrono>
#include <boost/program_options.hpp>

#include "../../src/graph/GraphIO.h"
#include "../../src/options.h"

#include "../../src/lower_bound/lower_bound_utils.h"
#include "../../src/version.h"


void
write_output_file(const std::string &path, const Configuration &config, const Instance &instance, const std::vector<Cost> &lower_bound_values,
                  const std::vector<double> &initialization_times, const std::vector<double> &result_times,
                  const std::vector<double> &complete_times, bool print_stdout = false) {
    using namespace YAML;

    Emitter out;
    out << BeginDoc << BeginMap;
    out << Key << "forbidden_subgraphs" << Value << config.forbidden_subgraphs;
    out << Key << "lower_bound_name" << Value << config.lower_bound;
    out << Key << "instance" << Value << instance.name;
    out << Key << "multiplier" << Value << instance.multiplier;
    out << Key << "permutation" << Value << instance.permutation;

    out << Key << "values" << Value << Flow << lower_bound_values;
    out << Key << "initialization_times" << Value << Flow << initialization_times << Comment("ns");
    out << Key << "result_times" << Value << Flow << result_times << Comment("ns");
    out << Key << "complete_times" << Value << Flow << complete_times << Comment("ns");

    out << Key << "commit_hash" << Value << GIT_COMMIT_HASH;
    out << Key << "timelimit" << Value << config.timelimit << Comment("s");
    out << EndMap << EndDoc;

    if (print_stdout)
        std::cout << out.c_str() << "\n";

    if (!path.empty()) {
        std::ofstream file(path);

        if (!file)
            throw std::runtime_error("could not open output file");

        file << out.c_str();
        file.close();
    }
}

int main(int argc, char *argv[]) {
    namespace po = boost::program_options;
    using Options::LB;

    // std::vector<Options::LB> lower_bounds = {LB::Greedy, LB::No, LB::LocalSearch, LB::LinearProgram};
    std::vector<std::string> inputs = {
            "../data/bio/bio-nr-3-size-16.graph",
            "../data/bio/bio-nr-4-size-39.graph",
            "../data/bio/bio-nr-11-size-22.graph",
            "../data/misc/karate.graph"
    };

    constexpr auto max_k = std::numeric_limits<Cost>::max();

    Configuration config(Options::FSG::C4P4, 100, Options::SolverType::FPT, Options::Selector::FirstFound, LB::Trivial);
    config.input_path = inputs[0];
    config.lower_bound = LB::SortedGreedy;
    config.permutation = 0;
    config.timelimit = -1;
    size_t iterations = 2;


    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("F", po::value<Options::FSG>(&config.forbidden_subgraphs)->default_value(config.forbidden_subgraphs))
            ("lower-bound", po::value<Options::LB>(&config.lower_bound)->default_value(config.lower_bound))
            ("input", po::value<std::string>(&config.input_path)->default_value(config.input_path),
             "path to input instance")
            ("output", po::value<std::string>(&config.output_path)->default_value(config.output_path), "output file")
            ("iterations", po::value<size_t>(&iterations)->default_value(iterations), "number of repetitions")
            ("permutation", po::value<int>(&config.permutation)->default_value(config.permutation),
             "permutation of input instance")
            ("multiplier", po::value<double>(&config.multiplier)->default_value(config.multiplier),
             "multiplier for discretization of input instance")
            ("timelimit", po::value<int>(&config.timelimit)->default_value(config.timelimit));

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }




    auto instance = GraphIO::read_instance(config.input_path, config.multiplier, config.permutation);
    VertexPairMap<bool> marked(instance.graph.size());


    write_output_file(config.output_path, config, instance, {}, {}, {}, {});



    constexpr auto FSG = Options::FSG::C4P4;

    if (config.forbidden_subgraphs != FSG) {
        throw std::runtime_error("Only C4P4 are currently supported.");
    }

    SubgraphStats<FSG> subgraph_stats(instance, marked);
    subgraph_stats.initialize(max_k);
    auto lb = lower_bound::make<FSG>(config.lower_bound, instance, marked, subgraph_stats, config);


    std::vector<double> initialization_times, result_times, complete_times;
    std::vector<Cost> values;
    for (size_t it = 0; it < iterations; ++it) {
        using namespace std::chrono;

        auto t1 = steady_clock::now();
        lb->initialize(max_k);
        auto t2 = steady_clock::now();
        auto value = lb->calculate_lower_bound(max_k);
        auto t3 = steady_clock::now();

        initialization_times.push_back(duration_cast<nanoseconds>(t2 - t1).count());
        result_times.push_back(duration_cast<nanoseconds>(t3 - t2).count());
        complete_times.push_back(duration_cast<nanoseconds>(t3 - t1).count());
        values.push_back(value);
    }

    write_output_file(config.output_path, config, instance, values, initialization_times, result_times, complete_times, true);

    return 0;
}

