#include <iostream>
#include <chrono>

#include <boost/program_options.hpp>

#include "graph/Graph.h"
#include "graph/GraphIO.h"

#include "Statistics.h"
#include "Editor.h"
#include "Permutation.h"
#include "Solution.h"

#include "finder/Center.h"
#include "finder/CenterC4P4.h"

#include "graph/synthetic_graphs.h"
#include "version.h"

#include "lower_bound/LinearProgramLowerBound.h"

#include "graph/algorithms.h"


void read_args(int argc, char *argv[], int& seed, std::string &input, double &multiplier, Cost &k_max, Options::FSG &fsg, Options::Selector &selector, Options::LB &lower_bound) {
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("input", po::value<std::string>(&input)->default_value(input), "path to input instance")
            ("multiplier", po::value<double>(&multiplier)->default_value(multiplier), "multiplier for discretization of input instances")
            ("seed", po::value<int>(&seed)->default_value(seed), "seed for instance permutation")
            ("k", po::value<Cost>(&k_max)->default_value(k_max), "maximum editing cost")
            ("F", po::value<Options::FSG>(&fsg)->default_value(fsg), "forbidden subgraphs")
            ("selector", po::value<Options::Selector>(&selector)->default_value(selector))
            ("lower_bound", po::value<Options::LB>(&lower_bound)->default_value(lower_bound))
            ;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        abort();
    }
}


int main(int argc, char *argv[]) {

    const std::vector<std::string> paths {
        "./data/bio/bio-nr-3-size-16.metis",
        "./data/bio/bio-nr-4-size-39.metis",
        "./data/bio/bio-nr-11-size-22.metis",
        "./data/karate.graph"};

     /*

    int seed = 0;
    std::string input = paths[2];
    double multiplier = 100;
    Cost k_max = 0;
    auto forbidden_type = Options::FSG::P4C4;
    auto selector = Options::Selector::FirstEditable;
    auto lower_bound = Options::LB::Greedy;

    read_args(argc, argv, seed, input, multiplier, k_max, forbidden_type, selector, lower_bound);


    auto instance = GraphIO::read_instance(input, multiplier);

    Permutation P(instance.graph.size(), seed);
    Permutation P_r = P.reverse();
    instance = P[instance];

    // Configuration config(0, Options::Selector::FirstEditable, Options::FSG::P4C4, Options::LB::No, "", "");
    // config.read_input(argc, argv);


    std::vector<Solution> solutions;

    Editor editor(instance, selector, forbidden_type, lower_bound);

    auto solution_cb = [&](const std::vector<VertexPair> &edits) {
        solutions.emplace_back(P_r[instance], P_r[edits]);
    };
    //auto pruning_cb = [](Cost k, Cost lb) { std::cout << "pruned: k=" << k << ", lb=" << lb << ", eps=" << lb - k << "\n"; };
    auto pruning_cb_2 = [](Cost, Cost) {};


    bool solved = editor.edit(k_max, solution_cb, pruning_cb_2);


    std::cout << (solved ? "instance solved" : "instance not solved") << "\n";

    std::sort(solutions.begin(), solutions.end(),
              [](const auto &lhs, const auto &rhs) { return lhs.cost > rhs.cost; });


    for (const auto& solution : solutions) {
        assert(solution.is_valid(P_r[instance], forbidden_type));
    }



    auto now = std::chrono::system_clock::now();
    auto itt = std::chrono::system_clock::to_time_t(now);

    std::ostringstream time_ss;
    time_ss << std::put_time(std::localtime(&itt), "%FT%TZ%z");

    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << "instance_name"
        << YAML::Value << instance.name;
    out << YAML::Key << "seed"
        << YAML::Value << seed;
    out << YAML::Key << "commit_hash"
        << YAML::Value << GIT_COMMIT_HASH;
    out << YAML::Key << "time"
        << YAML::Value << time_ss.str();
    out << YAML::Key << "algorithm"
        << YAML::Value << "fpt";
    out << YAML::Key << "forbidden_subgraphs"
        << YAML::Key << forbidden_type;
    out << YAML::Key << "config"
        << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "selector"
        << YAML::Value << selector;
    out << YAML::Key << "lower_bound"
        << YAML::Key << lower_bound;
    out << YAML::EndMap;
    out << YAML::Key << "solutions"
        << YAML::Value << solutions;
    out << YAML::Key << "statistics"
        << YAML::Value << editor.stats();
    out << YAML::EndMap;

    std::cout << "\n" << out.c_str() << "\n";

    auto print_subgraph = [](Subgraph &&subgraph) {
        std::cout << subgraph << "\n";
        return false;
    };

    auto adj_subgraphs = [](const Graph &G, VertexPair uv) {
        std::vector<Subgraph> subgraphs;
        detail::Center<4>(G).find_near(uv, [&](Subgraph &&subgraph) { subgraphs.push_back(std::move(subgraph)); return false; });
        return subgraphs;
    };

    Graph G(9);
    G.setEdges({{0, 1}, {1, 2}, {2, 3}, {1, 4}, {2, 5}, {4, 5}, {5, 6}, {6, 7}, {6, 8}});
    // 3-2-5-6-7
    //   | | |
    // 0-1-4 8

    for (VertexPair uv : G.vertexPairs()) {
        detail::Center<4>(G).find_near(uv, print_subgraph);
    }
    */


    // create_synthetic_cluster_graph(100, 0);



    VertexPairMap<Cost> map(10);

    std::vector<Solution> solutions = {{map, {{0, 1}}}, {map, {{0, 1}, {2, 3}}}, {map, {{0, 2}, {2, 4}}}};

    for (const auto &solution : solutions)
        std::cout << solution << " ";

    std::cout << "\n";

    Solution::filter_inclusion_minimal(solutions);
    for (const auto &solution : solutions)
        std::cout << solution << " ";


    auto instance = GraphIO::read_instance(paths[0], 100);
    VertexPairMap<bool> marked(instance.graph.size());
    std::shared_ptr<FinderI> finder = std::make_shared<Finder::CenterC4P4>(instance.graph);
    auto lb = LinearProgramLowerBound(instance, marked, finder);
    lb.initialize(1000);
    std::cout << lb.result(1000);

    return 0;
}