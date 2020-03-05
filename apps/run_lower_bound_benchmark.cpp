//
// Created by jonas on 04.09.19.
//


#include <iostream>
#include <chrono>
#include <boost/program_options.hpp>

#include "../src/graph/GraphIO.h"
#include "../src/options.h"

#include "../src/lower_bound/lower_bound_utils.h"
#include "../src/finder/finder_utils.h"


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;
    using Options::LB;

    // std::vector<Options::LB> lower_bounds = {LB::Greedy, LB::No, LB::LocalSearch, LB::LinearProgram};
    LB lower_bound = LB::LSSWZ_MWIS_Solver;
    std::vector<std::string> inputs = {
        "../data/bio/bio-nr-3-size-16.graph",
        "../data/bio/bio-nr-4-size-39.graph",
        "../data/bio/bio-nr-11-size-22.graph",
        "../data/misc/karate.graph"
    };
    std::string output;
    std::vector<int> seeds = {0, 1, 2, 3};
    size_t iterations = 10;
    double multiplier = 100;


    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("lower_bound", po::value<Options::LB>(&lower_bound)->default_value(lower_bound))
            ("inputs", po::value<std::vector<std::string>>()->multitoken(), "path to input instances")
            ("output", po::value<std::string>(&output)->default_value(output), "output file")
            ("iterations", po::value<size_t>(&iterations)->default_value(iterations), "number of repetitions")
            ("seeds", po::value<std::vector<int>>()->multitoken(), "seeds for permutation of instances")
            ("multiplier", po::value<double>(&multiplier)->default_value(multiplier), "multiplier for discretization of input instances")
            ;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }


    if (vm.count("inputs"))
        inputs = vm["inputs"].as<std::vector<std::string>>();

    if (vm.count("seeds"))
        seeds = vm["seeds"].as<std::vector<int>>();

    std::ofstream file;
    if (!output.empty()) {
        file = std::ofstream(output);
    }

    constexpr Cost max_k = std::numeric_limits<Cost>::max();


    Configuration config(Options::FSG::C4P4, multiplier, Options::SolverType::FPT, Options::Selector::FirstFound, lower_bound);


    for (const auto &input : inputs) {
        config.input_path = input;
    
        auto orig_instance = GraphIO::read_instance(input, config.multiplier);
        VertexPairMap<bool> marked(orig_instance.graph.size());
    
        
        for (auto seed : seeds) {
            config.seed = seed;

            Permutation P(orig_instance.graph.size(), config.seed);
            auto instance = P[orig_instance];

            std::shared_ptr<FinderI> finder = Finder::make(config.forbidden_subgraphs);
            SubgraphStats subgraph_stats(finder, instance, marked);
            subgraph_stats.initialize(max_k);
            auto lb = lower_bound::make(config.lower_bound, finder, instance, marked, subgraph_stats, config);


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
            
            using namespace YAML;
            Emitter out;
            out << BeginDoc << BeginMap;
            out << Key << "lower_bound_name" << Value << lower_bound;
            out << Key << "instance" << Value << instance.name;
            out << Key << "seed" << Value << seed;
            out << Key << "values" << Value << Flow << values;
            out << Key << "initialization_times" << Value << Flow << initialization_times << Comment("ns");
            out << Key << "result_times" << Value << Flow << result_times << Comment("ns");
            out << Key << "complete_times" << Value << Flow << complete_times << Comment("ns");
            out << EndMap << EndDoc;

            if (output.empty()) {
                std::cout << out.c_str();
            } else {
                file << out.c_str();
            }
        }
    }

    return 0;
}

