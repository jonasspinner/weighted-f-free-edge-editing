//
// Created by jonas on 19.08.19.
//


#include <boost/program_options.hpp>

#include "../src/solvers/FPTSolver.h"
#include "../src/graph/GraphIO.h"


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    const std::vector<std::string> paths {
            "../data/bio/bio-nr-3-size-16.metis",
            "../data/bio/bio-nr-11-size-22.metis",
            "../data/bio/bio-nr-4-size-39.metis",
            "../data/bio/bio-nr-915-size-68.metis",
            "../data/bio/bio-nr-936-size-71.metis",
            "../data/bio/bio-nr-970-size-96.metis",
            "../data/bio/bio-nr-277-size-222.metis",
            "../data/karate.graph",
            "../data/lesmis.graph",
            "../data/dolphins.graph",
            "../data/grass_web.metis.graph"};


    std::string input = "../data/bio/bio-nr-1020-size-44.metis";
    double multiplier = 100;
    Cost k_max = 6730;

    Configuration config(k_max, Options::Selector::MostMarked, Options::FSG::P4C4, Options::LB::LocalSearch,
            input, "", multiplier, Options::SolverType::FPT, 0, 0, false);

    auto options = config.options({Options::SolverType::FPT});

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(options).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << options << "\n";
        return 1;
    }


    auto instance = GraphIO::read_instance(config.input_path, config.multiplier);

    FPTSolver solver(config);
    auto result = solver.solve(instance);

    for (const auto &solution : result.solutions)
        std::cout << solution << "\n";

    return 0;
}