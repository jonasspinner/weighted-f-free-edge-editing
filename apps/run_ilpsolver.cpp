//
// Created by jonas on 17.08.19.
//


#include <vector>
#include <string>
#include <chrono>

#include <boost/program_options.hpp>

#include "../src/solvers/ILPSolver.h"
#include "../src/graph/GraphIO.h"


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;
    using namespace std::chrono;

    const std::vector<std::string> paths {
            "../data/bio/bio-nr-3-size-16.graph",
            "../data/bio/bio-nr-11-size-22.graph",
            "../data/bio/bio-nr-4-size-39.graph",
            "../data/bio/bio-nr-915-size-68.graph",
            "../data/bio/bio-nr-936-size-71.graph",
            "../data/bio/bio-nr-970-size-96.graph",
            "../data/bio/bio-nr-277-size-222.graph",
            "../data/misc/karate.graph",
            "../data/misc/lesmis.graph",
            "../data/misc/dolphins.graph",
            "../data/misc/grass_web.graph"};

    std::string input = paths[2];
    double multiplier = 100;

    Configuration config(Options::FSG::C4P4, input, multiplier, Options::SolverType::ILP, 0, Options::Selector::MostMarkedPairs, Options::LB::LocalSearch);

    auto options = config.options({Options::SolverType::ILP});

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(options).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << options << "\n";
        return 1;
    }

    std::ofstream outputFile;
    if (!config.output_path.empty())
        outputFile.open(config.output_path);


    auto instance = GraphIO::read_instance(config);

    ILPSolver solver(config);


    auto t1 = steady_clock::now();

    auto result = solver.solve(instance);

    auto t2 = steady_clock::now();
    auto solve_time = duration_cast<nanoseconds>(t2 - t1).count();


    auto solutions = result.solutions;
    std::sort(solutions.begin(), solutions.end());

    Cost solution_cost = -1;
    if (!solutions.empty())
        solution_cost = solutions[0].cost;


    using namespace YAML;
    Emitter out;
    out << BeginDoc << BeginMap;
    out << Key << "config" << Value << BeginMap;
    out << Key << "extended_constraints" << Value << config.extended_constraints;
    out << Key << "sparse_constraints" << Value << config.sparse_constraints;
    out << Key << "num_threads" << Value << config.num_threads;
    out << Key << "timelimit"<< Value << config.timelimit;
    out << Key << "lazy" << Value << config.lazy;
    out << EndMap;
    out << Key << "instance" << Value << instance;
    out << Key << "forbidden_subgraphs" << Value << config.forbidden_subgraphs;
    out << Key << "solutions" << Value << solutions;
    out << Key << "solution_cost" << Value << solution_cost;
    out << Key << "time" << Value << solve_time << Comment("ns");
    out << EndMap << EndDoc;

    std::cout << out.c_str() << "\n";
    if (outputFile)
        outputFile << out.c_str();

    return 0;
}