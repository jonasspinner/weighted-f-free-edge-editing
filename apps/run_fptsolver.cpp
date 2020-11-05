#include <chrono>
#include <boost/program_options.hpp>

#include "../src/solvers/FPTSolver.h"
#include "../src/graph/GraphIO.h"
#include "../src/version.h"


void write_output_file(const std::string &path, const Configuration &config, const Instance &instance,
                       const std::vector<Solution> &solutions, Cost solution_cost, std::chrono::nanoseconds solve_time,
                       const std::vector<FPTSolver::Stat> &stats, bool print_stdout = false) {
    using namespace YAML;

    std::vector<Cost> stats_k;
    std::vector<int> stats_calls;
    std::vector<double> stats_time;

    for (const auto &[k, c, t] : stats) {
        stats_k.push_back(k);
        stats_calls.push_back(c);
        stats_time.push_back(t);
    }

    Emitter out;
    out << BeginDoc << BeginMap;
    out << Key << "config" << Value << BeginMap;
    out << Key << "k_max" << Value << config.k_max;
    out << Key << "selector" << Value << config.selector;
    out << Key << "lower_bound" << Value << config.lower_bound;
    out << Key << "find_all_solutions" << Value << config.find_all_solutions;
    out << Key << "pre_mark_vertex_pairs" << Value << config.pre_mark_vertex_pairs;
    out << Key << "search_strategy" << Value << config.search_strategy;
    out << EndMap;
    out << Key << "stats" << Value << BeginMap;
    out << Key << "k" << Value << Flow << stats_k;
    out << Key << "calls" << Value << Flow << stats_calls;
    out << Key << "time" << Value << Flow << stats_time << Comment("ns");
    out << EndMap;
    out << Key << "commit_hash" << Value << GIT_COMMIT_HASH;
    out << Key << "instance" << Value << instance;
    out << Key << "forbidden_subgraphs" << Value << config.forbidden_subgraphs;
    out << Key << "solutions" << Value << solutions;
    out << Key << "solution_cost" << Value << solution_cost;
    out << Key << "time" << Value << solve_time.count() << Comment("ns");
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


    double multiplier = 100;

    Configuration config(Options::FSG::C4P4, multiplier, Options::SolverType::FPT, Options::Selector::MostAdjacentSubgraphs, Options::LB::Trivial);

    config.find_all_solutions = true;
    config.pre_mark_vertex_pairs = false;

    config.verbosity = 2;

    config.selector = Options::Selector::MostAdjacentSubgraphs;
    config.lower_bound = Options::LB::WeightedPackingLocalSearch;
    config.search_strategy = Options::FPTSearchStrategy::Fixed;
    //config.input_path = paths[0]; config.k_max = 0;
    //config.input_path = paths[0]; config.k_max = 426;
    //config.input_path = paths[1]; config.k_max = 1253;
    //config.input_path = "../data/bio/bio-nr-1020-size-44.graph"; config.k_max = 6730;
    // config.input_path = "../data/bio/bio-nr-1585-size-88.graph"; config.k_max = 56236; // 863 s
    // config.find_all_solutions = false;
    //config.input_path = "../data/bio/bio-nr-915-size-68.graph"; config.k_max = 14486 - 0; //14257;
    // config.input_path = "../data/bio/bio-nr-613-size-52.graph"; config.k_max = 22567;
    config.input_path = "../data/bio/bio-nr-782-size-17.graph"; config.k_max = 3833;
    config.input_path = "../data/bio/bio-nr-167-size-64.graph"; config.k_max = 13972;
    config.input_path = "../data/bio/bio-nr-1793-size-50.graph"; config.k_max = 18141;
    config.input_path = "../data/bio/bio-nr-32-size-50.graph"; config.k_max = 21184;
    //config.input_path = "../data/bio/bio-nr-1936-size-68.graph"; config.k_max = 43161;

    config.input_path = "../data/bio/bio-nr-1020-size-44.graph"; config.k_max = 6730;

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


    auto instance = GraphIO::read_instance(config);

    write_output_file(config.output_path, config, instance, {}, -1, nanoseconds(-1), {});

    if (config.search_strategy == Options::FPTSearchStrategy::Fixed && config.k_max == -1)
        return 0;

    FPTSolver solver(config);


    auto t1 = steady_clock::now();

    auto result = solver.solve(instance);

    auto t2 = steady_clock::now();
    nanoseconds solve_time(t2 - t1);

    auto solutions = result.solutions;
    std::sort(solutions.begin(), solutions.end());

    Cost solution_cost = -1;
    if (!solutions.empty())
        solution_cost = solutions[0].cost;
    

    write_output_file(config.output_path, config, instance, solutions, solution_cost, solve_time, solver.stats(), true);

    return 0;
}