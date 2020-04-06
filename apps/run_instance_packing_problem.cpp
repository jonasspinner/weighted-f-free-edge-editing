//
// Created by jonas on 06.04.20.
//


#include <boost/program_options.hpp>
#include <iostream>

#include "../src/graph/GraphIO.h"
#include "../src/finder/finder_utils.h"


std::vector<Subgraph> get_all_forbidden_subgraphs(const Graph &graph, Options::FSG fsg) {
    auto finder = Finder::make(fsg);
    std::vector<Subgraph> result;
    finder->find(graph, [&](Subgraph &&subgraph) {
        result.push_back(std::move(subgraph));
        return false;
    });
    return result;
}


std::pair<std::vector<VertexPair>, std::vector<Cost>>
get_covered_edges_and_costs(const std::vector<Subgraph> &subgraphs, const VertexPairMap<Cost> &costs) {
    Graph covered(costs.size());
    for (const auto &subgraph : subgraphs) {
        for (auto uv : subgraph.vertexPairs()) {
            covered.setEdge(uv);
        }
    }

    std::vector<VertexPair> edges;
    std::vector<Cost> edges_costs;
    for (auto uv : covered.edges()) {
        edges.push_back(uv);
        edges_costs.push_back(costs[uv]);
    }
    return {edges, edges_costs};
}


int main(int argc, char *argv[]) {
    namespace po = boost::program_options;

    try {

        std::string input_path;
        double multiplier = 1;
        int permutation = 0;
        std::string output_path;
        auto fsg = Options::FSG::C4P4;

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
                ("fsg", po::value<Options::FSG>(&fsg)->default_value(fsg));

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

        std::ostream *os = &std::cout;
        std::ofstream output_file;

        if (!output_path.empty()) {
            output_file.open(output_path);
            if (!output_file)
                throw std::runtime_error("could not open output_file");
            os = &output_file;
        }

        auto instance = GraphIO::read_instance(input_path, multiplier, permutation);

        auto forbidden_subgraphs = get_all_forbidden_subgraphs(instance.graph, fsg);
        auto[covered_edges, covered_edges_costs] = get_covered_edges_and_costs(forbidden_subgraphs, instance.costs);

        using namespace YAML;

        Emitter out;
        out << BeginDoc << BeginMap;
        out << Key << "forbidden_subgraphs" << Value << Flow << forbidden_subgraphs << Comment("conflict set");
        out << Key << "covered_edges" << Value << Flow << covered_edges;
        out << Key << "covered_edges_costs" << Value << Flow << covered_edges_costs;
        out << EndMap << EndDoc;

        *os << out.c_str() << "\n";


    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "unknown exception\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
