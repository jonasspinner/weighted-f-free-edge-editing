//
// Created by jonas on 06.04.20.
//


#include <boost/program_options.hpp>
#include <iostream>

#include "../src/graph/GraphIO.h"


template <typename T>
std::vector<std::vector<T>> to_matrix(const VertexPairMap<T> &map) {
    auto n = map.size();
    auto out = std::vector<std::vector<T>>(n, std::vector<T>(n));

    for (auto uv : Graph::VertexPairs(n)) {
        out[uv.u][uv.v] = map[uv];
        out[uv.v][uv.u] = map[uv];
    }

    return out;
}


std::vector<std::vector<int>> to_matrix(const Graph &graph) {
    auto n = graph.size();
    auto out = std::vector<std::vector<int>>(n, std::vector<int>(n));

    for (auto uv : graph.edges()) {
        out[uv.u][uv.v] = 1;
        out[uv.v][uv.u] = 1;
    }

    return out;
}


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    try {

        std::string input_path;
        double multiplier = 1;
        int permutation = 0;
        std::string output_path = "-";

        po::options_description desc("Allowed options");
        desc.add_options()
                ("help", "produce help message")
                ("input", po::value<std::string>(&input_path)->required(), "path to input instance")
                ("multiplier", po::value<double>(&multiplier)->default_value(multiplier))
                ("permutation", po::value<int>(&permutation)->default_value(permutation))
                ("output", po::value<std::string>(&output_path)->default_value(output_path))
                ;

        po::positional_options_description pd;
        pd.add("input", 1);
        pd.add("output", 2);


        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(pd).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 1;
        }



        auto instance = GraphIO::read_instance(input_path, multiplier, permutation);

        using namespace YAML;

        Emitter out;
        out << BeginDoc << BeginMap;
        out << Key << "instance" << Value << input_path;
        out << Key << "multiplier" << Value << multiplier;
        out << Key << "permutation" << Value << permutation;
        out << Key << "adj" << Value << Flow << to_matrix(instance.graph);
        out << Key << "costs" << Value << Flow << to_matrix(instance.costs);
        out << EndMap << EndDoc;

        if (!output_path.empty()) {
            if (output_path == "-") {
                std::cout << out.c_str();
            } else {
                GraphIO::write_instance(output_path, instance.graph, instance.costs);
            }
        }
    } catch (const std::exception &e) {
        std::cerr << e.what();
        return 1;
    } catch (...) {
        std::cerr << "unknown exception\n";
        return 1;
    }

    return 0;
}