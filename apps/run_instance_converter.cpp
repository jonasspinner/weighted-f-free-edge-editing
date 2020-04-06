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


void to_yaml(const Instance &instance, const std::string& input_path, double multiplier, int permutation, std::ostream& os) {
    using namespace YAML;

    Emitter out;
    out << BeginDoc << BeginMap;
    out << Key << "instance" << Value << input_path;
    out << Key << "multiplier" << Value << multiplier;
    out << Key << "permutation" << Value << permutation;
    out << Key << "adj" << Value << Flow << to_matrix(instance.graph);
    out << Key << "costs" << Value << Flow << to_matrix(instance.costs);
    out << EndMap << EndDoc;

    os << out.c_str();
}


// Note: Other formats may include metis.
enum class OutputFormat {
    YAML
};

std::istream& operator>>(std::istream& in, OutputFormat& selector) {
    std::string token;
    in >> token;
    if (token == "yaml")
        selector = OutputFormat::YAML;
    else
        in.setstate(std::ios_base::failbit);
    return in;
}

std::ostream &operator<<(std::ostream &os, OutputFormat selector) {
    switch (selector) {
        case OutputFormat::YAML:
            return os << "yaml";
        default:
            return os;
    }
}


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;

    try {

        std::string input_path;
        double multiplier = 1;
        int permutation = 0;
        std::string output_path;

        auto output_type = OutputFormat::YAML;

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "Produce help message.")
            ("input", po::value<std::string>(&input_path)->required(),
                "Path to input instance.\nMay be provided as named option or as the first positional "
                "parameter.")
            ("multiplier", po::value<double>(&multiplier)->default_value(multiplier),
                "Multiplier for multiplying editing costs before discretization.")
            ("permutation", po::value<int>(&permutation)->default_value(permutation),
                "Permutation seed for the instance. The value 0 represents the identity.")
            ("output", po::value<std::string>(&output_path)->default_value(output_path),
                "Path to output file. If no path is specified the output is written to stdout.\n"
                "May be provided as named option or as the second positional parameter.")
            ("format", po::value<OutputFormat>(&output_type)->default_value(output_type),
                "Output format. Currently only 'yaml' is possible.")
            ;

        po::positional_options_description pd;
        pd.add("input", 1);
        pd.add("output", 2);


        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(pd).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << "Usage: " << argv[0] << " [options...] input [output]\n";
            std::cout << desc << "\n";
            return EXIT_FAILURE;
        }

        std::ostream* os = &std::cout;
        std::ofstream output_file;

        if (!output_path.empty()) {
            output_file.open(output_path);
            if (!output_file)
                throw std::runtime_error("could not open output_file");
            os = &output_file;
        }

        auto instance = GraphIO::read_instance(input_path, multiplier, permutation);

        switch (output_type) {
            case OutputFormat::YAML:
                to_yaml(instance, input_path, multiplier, permutation, *os);
                break;
            default:
                throw std::runtime_error("unknown value for output_type");
        }

    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "unknown exception\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}