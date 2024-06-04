//
// Created by jonas on 12.09.19.
//


#include <chrono>
#include "../src/Configuration.h"
#include "../src/editor/Editor.h"
#include "../src/graph/GraphIO.h"
#include "../src/Solution.h"

int main(int argc, char* argv[]) {
    namespace po = boost::program_options;
    using namespace std::chrono;

    const std::vector<std::string> paths {
            "../data/test/bio-nr-3-size-16.graph",
            "../data/bio/bio-nr-3-size-16.graph",
            "../data/bio/bio-nr-4-size-39.graph",
            "../data/bio/bio-nr-11-size-22.graph",
            "../data/bio/bio-nr-277-size-222.graph",
            "../data/misc/karate.graph",
            "../data/misc/lesmis.graph",
            "../data/misc/dolphins.graph",
            "../data/misc/grass_web.graph"};

    std::string input = paths[0];
    double multiplier = 100;
    Cost k_max = 10780;
    int num_steps = 10;
    int num_steps_warming_up = 10;
    int timeout = 100;

    po::options_description additional_options("Additional options");
    additional_options.add_options()
            ("help", "produce help message")
            ("input", po::value<std::string>(&input)->default_value(input), "path to input instance")
            ("multiplier", po::value<double>(&multiplier)->default_value(multiplier), "multiplier for discretization of input instances")
            ("k", po::value<Cost>(&k_max)->default_value(k_max), "maximum editing cost")
            ("steps", po::value<int>(&num_steps)->default_value(num_steps), "")
            ("warmup-steps", po::value<int>(&num_steps_warming_up)->default_value(num_steps_warming_up), "")
            ("timeout", po::value<int>(&timeout)->default_value(timeout))
        // ("output", po::value<std::string>(), "path to output file")
            ;

    po::variables_map vm;
    try {
        auto parsed_options = po::command_line_parser(argc, argv).options(additional_options).run();
        po::store(parsed_options, vm);
        po::notify(vm);
    } catch (const po::invalid_option_value &e) {
        std::cerr << "invalid options value " << e.get_option_name() << "\n";
        return 1;
    } catch (const po::unknown_option &e) {
        std::cerr << "unknown option " << e.get_option_name() << "\n";
        return 1;
    }

    if (vm.count("help")) {
        std::cout << additional_options << "\n";
        return 1;
    }

    Configuration config(Options::FSG::C4P4, multiplier, Options::SolverType::FPT, Options::Selector::FirstFound, Options::LB::LocalSearch);
    config.find_all_solutions = false;
    config.input_path = input;
    config.k_max = k_max;

    // std::cout << config << "\n";

    auto instance = GraphIO::read_instance(config);
    Editor editor(instance.graph.copy(), instance.costs, config);

    Cost k_init = editor.initial_lower_bound();
    //for (int i = 0; i < num_steps_warming_up; ++i)
    //    editor.edit(0, [&](const std::vector<VertexPair> &){}, [](Cost, Cost){});

    std::vector<Cost> ks;
    std::vector<double> times;
    std::vector<Solution> solutions;
    bool timed_out = false;
    bool solved = false;

    for (int i = 0; i <= num_steps; ++i) {
        Cost k = k_init + (i * (config.k_max - k_init)) / num_steps;

        solutions.clear();

        auto t1 = std::chrono::steady_clock::now();

        solved = editor.edit(k, [&](const std::vector<VertexPair> &edits){ solutions.emplace_back(instance, edits); }, [](Cost, Cost){});

        auto t2 = std::chrono::steady_clock::now();
        auto time = duration_cast<nanoseconds>(t2 - t1).count();
        std::cout << "edit(" << k << ") took " << time << " ns\n";
        ks.push_back(k); times.push_back(time);

        if (!solutions.empty()) {
            break;
        } else if (timeout >= 0 && time > duration_cast<nanoseconds>(seconds(timeout)).count()) {
            timed_out = true;
            break;
        }
    }

    using namespace YAML;
    Emitter out;
    out << BeginMap;
    out << Key << "instance" << Value << instance;
    out << Key << "forbidden_subgraphs" << Value << config.forbidden_subgraphs;
    out << Key << "k_max" << Value << config.k_max;
    out << Key << "k" << Value << Flow << ks;
    out << Key << "time" << Value << Flow << times << Comment("ns");
    out << Key << "solutions" << Value << Flow << solutions;
    out << Key << "status" << Value << (timed_out ? "Timeout" : solved ? "Solved" : "Unsolved");
    out << EndMap;

    std::cout << out.c_str();
}