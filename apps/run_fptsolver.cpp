//
// Created by jonas on 19.08.19.
//


#include <chrono>
#include <boost/program_options.hpp>

#include "../src/solvers/FPTSolver.h"
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


    std::string input = "../data/bio/bio-nr-1020-size-44.graph";
    double multiplier = 100;
    Cost k_max = 6730;

    Configuration config(Options::FSG::C4P4, input, multiplier, Options::SolverType::FPT, k_max, Options::Selector::MostMarkedPairs, Options::LB::LocalSearch);

    auto options = config.options({Options::SolverType::FPT});

    po::variables_map vm;
    try {
        auto parsed_options = po::command_line_parser(argc, argv).options(options).run();
        po::store(parsed_options, vm);
        po::notify(vm);
    } catch (const po::invalid_option_value &e) {
        std::cerr << "invalid options value for option " << e.get_option_name() << "\n";
        return 1;
    } catch (const po::unknown_option &e) {
        std::cerr << "unknown option " << e.get_option_name() << "\n";
        return 1;
    } catch (const po::error &e) {
        std::cerr << "an error occured parsing the program options\n";
    }


    if (vm.count("help")) {
        std::cout << options << "\n";
        return 1;
    }

    std::ofstream outputFile;
    if (!config.output_path.empty())
        outputFile.open(config.output_path);


    auto instance = GraphIO::read_instance(config);

    FPTSolver solver(config);


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
    out << Key << "k_max" << Value << config.k_max;
    out << Key << "selector" << Value << config.selector;
    out << Key << "lower_bound" << Value << config.lower_bound;
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