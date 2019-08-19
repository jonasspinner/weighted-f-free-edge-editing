//
// Created by jonas on 17.08.19.
//


#include <vector>
#include <string>

#include "../src/ILPSolver.h"
#include "../src/graph/GraphIO.h"


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    const std::vector<std::string> paths {
            "../data/cost_matrix_component_nr_3_size_16_cutoff_10.0.metis",
            "../data/cost_matrix_component_nr_4_size_39_cutoff_10.0.metis",
            "../data/cost_matrix_component_nr_11_size_22_cutoff_10.0.metis",
            "../data/cost_matrix_component_nr_277_size_222_cutoff_10.0.metis",
            "./data/karate.graph",
            "./data/lesmis.graph",
            "./data/dolphins.graph",
            "./data/grass_web.metis.graph"};


    std::string input = paths[0];
    double multiplier = 100;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("input", po::value<std::string>(&input)->default_value(input), "path to input instance")
            ("multiplier", po::value<double>(&multiplier)->default_value(multiplier), "multiplier for discretization of input instances")
            ("output", po::value<std::string>(), "path to output file")
            ;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }


    auto instance = GraphIO::read_graph(input, multiplier);

    ILPSolver solver(Options::FSG::P4C4);
    auto solutions = solver.solve(instance);

    
    for (const auto &solution : solutions)
        std::cout << solution << "\n";

    return 0;
}