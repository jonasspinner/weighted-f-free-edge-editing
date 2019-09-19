//
// Created by jonas on 30.08.19.
//


#include <boost/program_options.hpp>
#include <chrono>

#include "../src/interfaces/FinderI.h"
#include "../src/finder/Center.h"
#include "../src/finder/CenterC4P4.h"
#include "../src/finder/CenterP3.h"
#include "../src/finder/NaiveC4P4.h"
#include "../src/finder/NaiveP3.h"
#include "../src/graph/GraphIO.h"
#include "../src/Permutation.h"


std::unique_ptr<FinderI> make_finder(const std::string &name, const Graph &graph) {
    if (name == "CenterRecC4P4") {
        return std::make_unique<Finder::CenterRecC4P4>(graph);
    } else if (name == "CenterC4P4") {
        return std::make_unique<Finder::CenterC4P4>(graph);
    } else if (name == "NaiveC4P4") {
        return std::make_unique<Finder::NaiveC4P4>(graph);
    } else if (name == "CenterRecP3") {
        return std::make_unique<Finder::CenterRecP3>(graph);
    } else if (name == "CenterP3") {
        return std::make_unique<Finder::CenterP3>(graph);
    } else if (name == "NaiveP3") {
        return std::make_unique<Finder::NaiveP3>(graph);
    } else {
        return nullptr;
    }
}


double mean(const std::vector<double> &X) {
    double sum = 0;
    for (auto x : X)
        sum += x;

    return sum / X.size();
}


double standard_deviation(const std::vector<double> &X) {
    auto m = mean(X);
    double sum = 0;
    for (auto x : X)
        sum += (x - m) * (x - m);

    return sqrt(sum / (X.size() - 1));
}


void find_all_subgraphs_benchmark(YAML::Emitter &out, FinderI &finder, size_t iterations, int seed) {
    std::vector<double> find_all_times(iterations);
    int count = 0;

    for (size_t it = 0; it < iterations; ++it) {
        count = 0;
        auto start = std::chrono::system_clock::now();

        finder.find([&](Subgraph &&) {
            ++count;
            return false;
        });

        auto end = std::chrono::system_clock::now();
        find_all_times[it] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    }

    out << YAML::BeginMap;
    out << YAML::Key << "name" << YAML::Value << "find all subgraphs";
    out << YAML::Key << "seed" << YAML::Value << seed;
    out << YAML::Key << "count" << YAML::Value << count;
    out << YAML::Key << "time" << YAML::Value << YAML::Flow << find_all_times << YAML::Comment("ns");
    out << YAML::Key << "time_mean" << YAML::Value << mean(find_all_times) << YAML::Comment("ns");
    out << YAML::Key << "time_std" << YAML::Value << standard_deviation(find_all_times)
        << YAML::Comment("ns");
    out << YAML::EndMap;
}

void find_one_subgraph_benchmark(YAML::Emitter &out, FinderI &finder, size_t iterations, int seed) {
    std::vector<double> find_one_times(iterations);
    int count = 0;

    for (size_t it = 0; it < iterations; ++it) {
        count = 0;
        auto start = std::chrono::system_clock::now();

        finder.find([&](Subgraph &&) {
            ++count;
            return true;
        });

        auto end = std::chrono::system_clock::now();
        find_one_times[it] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    }

    out << YAML::BeginMap;
    out << YAML::Key << "name" << YAML::Value << "find one subgraphs";
    out << YAML::Key << "seed" << YAML::Value << seed;
    out << YAML::Key << "count" << YAML::Value << count;
    out << YAML::Key << "time" << YAML::Value << YAML::Flow << find_one_times << YAML::Comment("ns");
    out << YAML::Key << "time_mean" << YAML::Value << mean(find_one_times) << YAML::Comment("ns");
    out << YAML::Key << "time_std" << YAML::Value << standard_deviation(find_one_times)
        << YAML::Comment("ns");
    out << YAML::EndMap;
}


void find_all_near_subgraphs(YAML::Emitter &out, FinderI &finder, const Graph &graph, size_t iterations, int seed) {
    std::vector<double> find_all_near_times(iterations);
    int count = 0;

    for (size_t it = 0; it < iterations; ++it) {
        count = 0;
        auto start = std::chrono::system_clock::now();

        for (VertexPair uv : graph.vertexPairs()) {
            finder.find_near(uv, [&](Subgraph &&) {
                ++count;
                return false;
            });
        }

        auto end = std::chrono::system_clock::now();
        find_all_near_times[it] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    }

    out << YAML::BeginMap;
    out << YAML::Key << "name" << YAML::Value << "find all near subgraphs";
    out << YAML::Key << "seed" << YAML::Value << seed;
    out << YAML::Key << "count" << YAML::Value << count;
    out << YAML::Key << "time" << YAML::Value << YAML::Flow << find_all_near_times << YAML::Comment("ns");
    out << YAML::Key << "time_mean" << YAML::Value << mean(find_all_near_times) << YAML::Comment("ns");
    out << YAML::Key << "time_std" << YAML::Value << standard_deviation(find_all_near_times)
        << YAML::Comment("ns");
    out << YAML::EndMap;
}


void find_one_near_subgraph(YAML::Emitter &out, FinderI &finder, const Graph &graph, size_t iterations, int seed) {
    std::vector<double> find_one_near_times(iterations);
    int count = 0;

    for (size_t it = 0; it < iterations; ++it) {
        count = 0;
        auto start = std::chrono::system_clock::now();

        for (VertexPair uv : graph.vertexPairs()) {
            finder.find_near(uv, [&](Subgraph &&) {
                ++count;
                return true;
            });
        }

        auto end = std::chrono::system_clock::now();
        find_one_near_times[it] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    }

    out << YAML::BeginMap;
    out << YAML::Key << "name" << YAML::Value << "find one near subgraphs";
    out << YAML::Key << "seed" << YAML::Value << seed;
    out << YAML::Key << "count" << YAML::Value << count;
    out << YAML::Key << "time" << YAML::Value << YAML::Flow << find_one_near_times << YAML::Comment("ns");
    out << YAML::Key << "time_mean" << YAML::Value << mean(find_one_near_times) << YAML::Comment("ns");
    out << YAML::Key << "time_std" << YAML::Value << standard_deviation(find_one_near_times)
        << YAML::Comment("ns");
    out << YAML::EndMap;
}


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    std::string finder_name = "CenterC4P4";
    std::vector<std::string> inputs = {"./data/bio/bio-nr-3-size-16.metis"};
    std::vector<int> seeds = {0, 1};
    size_t iterations = 100;
    double multiplier = 100;


    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("finder_name", po::value<std::string>(&finder_name)->default_value(finder_name))
            ("inputs", po::value<std::vector<std::string>>()->multitoken(), "path to input instances")
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



    for (const auto &input : inputs) {
        auto instance = GraphIO::read_instance(input, multiplier);

        YAML::Emitter out;
        out << YAML::BeginDoc << YAML::BeginMap;
        out << YAML::Key << "instance" << YAML::Value << instance.name;
        out << YAML::Key << "finder" << YAML::Value << finder_name;
        out << YAML::Key << "iterations" << YAML::Value << iterations;
        out << YAML::Key << "experiments" << YAML::Value << YAML::BeginSeq;

        for (auto seed : seeds) {
            Permutation P(instance.graph.size(), seed);
            auto graph = P[instance.graph];

            auto finder = make_finder(finder_name, graph);

            find_all_subgraphs_benchmark(out, *finder, iterations, seed);
            find_one_subgraph_benchmark(out, *finder, iterations, seed);

            // find_all_near_subgraphs(out, *finder, graph, iterations, seed);
            // find_one_near_subgraph(out, *finder, graph, iterations, seed);
        }

        out << YAML::EndSeq << YAML::EndMap << YAML::EndDoc;
        std::cout << out.c_str() << "\n";
    }

    return 0;
}