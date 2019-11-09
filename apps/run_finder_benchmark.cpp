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
#include "../src/finder/Endpoint.h"
#include "../src/graph/GraphIO.h"
#include "../src/version.h"
#include "../src/finder/Naive.h"
#include "../src/finder/OuterP3.h"


std::unique_ptr<FinderI> make_finder(const std::string &name, const Graph &graph) {
    if (name == "CenterRecC6P6") {
        return std::make_unique<Finder::CenterRecC6P6>(graph);
    } if (name == "CenterRecC5P5") {
        return std::make_unique<Finder::CenterRecC5P5>(graph);
    } else if (name == "CenterRecC4P4") {
        return std::make_unique<Finder::CenterRecC4P4>(graph);
    } else if (name == "CenterRecP6") {
        return std::make_unique<Finder::CenterRecP6>(graph);
    } else if (name == "CenterRecP5") {
        return std::make_unique<Finder::CenterRecP5>(graph);
    } else if (name == "CenterRecP4") {
        return std::make_unique<Finder::CenterRecP4>(graph);
    } else if (name == "CenterRecP3") {
        return std::make_unique<Finder::CenterRecP3>(graph);
    } else if (name == "EndpointRecC6P6") {
        return std::make_unique<Finder::EndpointRecC6P6>(graph);
    } else if (name == "EndpointRecC5P5") {
        return std::make_unique<Finder::EndpointRecC5P5>(graph);
    } else if (name == "EndpointRecC4P4") {
        return std::make_unique<Finder::EndpointRecC4P4 >(graph);
    } else if (name == "EndpointRecP6") {
        return std::make_unique<Finder::EndpointRecP6>(graph);
    } else if (name == "EndpointRecP5") {
        return std::make_unique<Finder::EndpointRecP5>(graph);
    } else if (name == "EndpointRecP4") {
        return std::make_unique<Finder::EndpointRecP4>(graph);
    } else if (name == "EndpointRecP3") {
        return std::make_unique<Finder::EndpointRecP3>(graph);
    } else if (name == "CenterC4P4") {
        return std::make_unique<Finder::CenterC4P4>(graph);
    } else if (name == "CenterP3") {
        return std::make_unique<Finder::CenterP3>(graph);
    } else if (name == "NaiveC4P4") {
        return std::make_unique<Finder::NaiveC4P4>(graph);
    } else if (name == "NaiveP3") {
        return std::make_unique<Finder::NaiveP3>(graph);
    } else if (name == "NaiveRecC6P6") {
        return std::make_unique<Finder::NaiveRecC6P6>(graph);
    } else if (name == "NaiveRecC5P5") {
        return std::make_unique<Finder::NaiveRecC5P5>(graph);
    } else if (name == "NaiveRecC4P4") {
        return std::make_unique<Finder::NaiveRecC4P4>(graph);
    } else if (name == "NaiveRecP6") {
        return std::make_unique<Finder::NaiveRecP6>(graph);
    } else if (name == "NaiveRecP5") {
        return std::make_unique<Finder::NaiveRecP4>(graph);
    } else if (name == "NaiveRecP4") {
        return std::make_unique<Finder::NaiveRecP3>(graph);
    } else if (name == "NaiveRecP3") {
        return std::make_unique<Finder::NaiveRecP3>(graph);
    } else if (name == "OuterP3") {
        return std::make_unique<Finder::OuterP3>(graph);
    } else {
        std::cerr << "name = " << name << "\n";
        throw std::runtime_error("Finder name not valid.");
    }
}

bool has_near(const std::string &name) {
    if (name == "CenterRecC5P5") {
        return false;
    } else if (name == "CenterRecC4P4") {
        return false;
    } else if (name == "CenterRecP3") {
        return false;
    } else if (name == "EndpointRecC5P5") {
        return false;
    } else if (name == "EndpointRecC4P4") {
        return false;
    } else if (name == "EndpointRecP3") {
        return false;
    } else if (name == "CenterC4P4") {
        return true;
    } else if (name == "CenterP3") {
        return true;
    } else if (name == "NaiveC4P4") {
        return true;
    } else if (name == "NaiveP3") {
        return true;
    }
    return false;
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

    return sqrt(sum / static_cast<double>(X.size() - 1));
}

void write_output_doc(std::ostream &file, const FinderI &finder, const std::string &type, const Instance &instance,
        int count, const std::vector<double> &times, bool print_stdout = false) {
    using namespace YAML;

    Emitter out;
    out << BeginDoc << BeginMap;
    out << Key << "type" << Value << type;
    out << Key << "finder" << Value << finder.name();
    out << Key << "forbidden_subgraphs" << Value << finder.forbidden_subgraphs();
    out << Key << "commit_hash" << Value << GIT_COMMIT_HASH;
    out << Key << "instance" << Value << instance;
    out << Key << "count" << Value << count;
    out << Key << "time" << Value << Flow << times << Comment("ns");
    out << Key << "time_mean" << Value << mean(times) << Comment("ns");
    out << Key << "time_std" << Value << standard_deviation(times)
        << Comment("ns");
    out << EndMap << EndDoc;

    if (print_stdout)
        std::cout << out.c_str() << "\n";

    if (file)
        file << out.c_str() << "\n";
}


std::pair<int, std::vector<double>> find_all_subgraphs_benchmark(FinderI &finder, size_t iterations) {
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

    return {count, find_all_times};
}

std::pair<int, std::vector<double>> find_one_subgraph_benchmark(FinderI &finder, size_t iterations) {
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

    return {count, find_one_times};
}


std::pair<int, std::vector<double>> find_all_near_subgraphs_benchmark(FinderI &finder, const Graph &graph, size_t iterations) {
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

    return {count, find_all_near_times};
}


std::pair<int, std::vector<double>> find_one_near_subgraph_benchmark(FinderI &finder, const Graph &graph, size_t iterations) {
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

    return {count, find_one_near_times};
}


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    std::string finder_name = "CenterC4P4";
    std::string input = "../data/bio/bio-nr-3-size-16.graph";
    std::string output_path;
    size_t iterations = 100;
    int permutation = 0;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("finder", po::value<std::string>(&finder_name)->default_value(finder_name))
            ("input", po::value<std::string>(&input)->default_value(input), "path to input instance")
            ("permutation", po::value<int>(&permutation)->default_value(permutation))
            ("output", po::value<std::string>(&output_path)->default_value(output_path))
            ("iterations", po::value<size_t>(&iterations)->default_value(iterations), "number of repetitions")
            ;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }


    // bool with_near = has_near(finder_name);
    bool with_near = false;

    std::ofstream file;
    if (!output_path.empty()) {
        file.open(output_path);

        if (!file)
            throw std::runtime_error("could not open output file");
    }

    auto original_instance = GraphIO::read_instance(input);
    Permutation P(original_instance.graph.size(), permutation);
    auto instance = P[original_instance];


    auto finder = make_finder(finder_name, instance.graph);

    write_output_doc(file, *finder, "find_all_subgraphs", instance, -1, {});
    write_output_doc(file, *finder, "find_one_subgraph", instance, -1, {});
    if (with_near) {
        write_output_doc(file, *finder, "find_all_near_subgraphs", instance, -1, {});
        write_output_doc(file, *finder, "find_one_near_subgraph", instance, -1, {});
    }
    file.close();


    int c1 = -1, c2 = -1, c3 = -1, c4 = -1;
    std::vector<double> t1, t2, t3, t4;

    std::tie(c1, t1) = find_all_subgraphs_benchmark(*finder, iterations);
    std::tie(c2, t2) = find_one_subgraph_benchmark(*finder, iterations);

    if (with_near) {
        std::tie(c3, t3) = find_all_near_subgraphs_benchmark(*finder, instance.graph, iterations);
        std::tie(c4, t4) = find_one_near_subgraph_benchmark(*finder, instance.graph, iterations);
    }


    if (!output_path.empty())
        file.open(output_path);

    write_output_doc(file, *finder, "find_all_subgraphs", instance, c1, t1, true);
    write_output_doc(file, *finder, "find_one_subgraph", instance, c2, t2, true);
    if (with_near) {
        write_output_doc(file, *finder, "find_all_near_subgraphs", instance, c3, t3, true);
        write_output_doc(file, *finder, "find_one_near_subgraph", instance, c4, t4, true);
    }
    return 0;
}