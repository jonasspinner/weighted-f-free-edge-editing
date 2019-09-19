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

    std::string input = paths[2];
    double multiplier = 100;

    Configuration config(0, Options::Selector::MostMarked, Options::FSG::P4C4, Options::LB::LocalSearch,
                         input, "", multiplier, Options::SolverType::ILP, 0, 0, false);

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

    auto instance = GraphIO::read_instance(config.input_path, config.multiplier);

    std::cout << "Solve " << instance.name << "\n";

    ILPSolver solver(config.forbidden_subgraphs);

    std::vector<Solution> solutions;
    std::vector<double> solve_times;
    for (int it = 0; it < 1; ++it) {
        auto t1 = steady_clock::now();
        auto result = solver.solve(instance);
        solutions = std::move(result.solutions);
        auto t2 = steady_clock::now();
        solve_times.push_back(duration_cast<microseconds>(t2 - t1).count());
    }


    std::sort(solutions.begin(), solutions.end());

    Cost solution_cost = -1;
    if (!solutions.empty())
        solution_cost = solutions[0].cost;

    YAML::Emitter out;
    out << YAML::BeginDoc << YAML::BeginMap;
    out << YAML::Key << "instance" << YAML::Value << instance;
    out << YAML::Key << "forbidden_subgraphs" << YAML::Value << config.forbidden_subgraphs;
    out << YAML::Key << "solutions" << YAML::Value << solutions;
    out << YAML::Key << "solution_cost" << YAML::Value << solution_cost;
    out << YAML::Key << "time" << YAML::Value << YAML::Flow << solve_times << YAML::Comment("ns");
    out << YAML::EndMap << YAML::EndDoc;

    std::cout << out.c_str() << "\n";
    if (outputFile)
        outputFile << out.c_str();

    return 0;
}