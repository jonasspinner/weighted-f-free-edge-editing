//
// Created by jonas on 19.08.19.
//


#include <boost/program_options.hpp>

#include "../src/FPTSolver.h"
#include "../src/graph/GraphIO.h"


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    const std::vector<std::string> paths {
            "./data/bio/bio-nr-3-size-16.metis",
            "./data/bio/bio-nr-4-size-39.metis",
            "./data/bio/bio-nr-11-size-22.metis",
            "./data/bio/bio-nr-277-size-222.metis",
            "./data/karate.graph",
            "./data/lesmis.graph",
            "./data/dolphins.graph",
            "./data/grass_web.metis.graph"};


    std::string input = paths[0];
    double multiplier = 100;
    Cost k_max = 600;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("input", po::value<std::string>(&input)->default_value(input), "path to input instance")
            ("multiplier", po::value<double>(&multiplier)->default_value(multiplier), "multiplier for discretization of input instances")
            ("k", po::value<Cost>(&k_max)->default_value(k_max), "maximum editing cost")
            // ("output", po::value<std::string>(), "path to output file")
            ;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }


    auto instance = GraphIO::read_graph(input, multiplier);

    FPTSolver solver(k_max, Options::FSG::P4C4);
    auto solutions = solver.solve(instance);

    for (const auto &solution : solutions)
        std::cout << solution << "\n";

    return 0;
}