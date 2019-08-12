#include <iostream>

#include "graph/Graph.h"
#include "graph/GraphIO.h"
#include "Statistics.h"
#include "Editor.h"

#include "finder/NaiveP3.h"
#include "finder/CenterC4P4.h"

#include "lower_bound/ExplicitMWISSolver.h"

#include "tests/FinderTests.h"
#include "tests/GraphTests.h"
#include "tests/SubgraphTests.h"
#include "tests/EditorTests.h"


class Solution {
public:
    Cost cost;
    std::vector<VertexPair> edits;

    // Solution(Cost cost, std::vector<VertexPair> edits) : cost(cost), edits(std::move(edits)) {}
    Solution(const VertexPairMap<Cost> &costs, std::vector<VertexPair> edits_) : cost(0), edits(std::move(edits_)) {
        for (VertexPair uv : edits) cost += costs[uv];
    }
    Solution(const Instance &instance, std::vector<VertexPair> edits_) : Solution(instance.costs, std::move(edits_)) {}

    friend std::ostream &operator<<(std::ostream &os, const Solution &solution) {
        os << solution.cost << " { ";
        for (auto uv : solution.edits) os << uv << " ";
        return os << "}";
    }
};

[[nodiscard]] bool is_solution_valid(const Graph &graph, const Solution &solution, Options::FSG forbidden) {
    Graph graph_(graph);
    std::unique_ptr<FinderI> finder;
    switch (forbidden) {
        case Options::FSG::P3:
            finder = std::make_unique<Finder::NaiveP3>(graph_);
            break;
        case Options::FSG::P4C4:
            finder = std::make_unique<Finder::CenterC4P4>(graph_);
            break;
    }

    for (VertexPair uv : solution.edits)
        graph_.toggle_edge(uv);

    bool found_forbidden_subgraph = finder->find([&](const Subgraph &) { return true; });

    return !found_forbidden_subgraph;
}


void search(const Instance &instance, const Configuration &config) {
    Cost k = 0;
    bool solved = false;
    do {
        Cost min_delta = std::numeric_limits<Cost>::max();
        std::vector<Cost> deltas;
        std::cout << "edit(" << std::setw(10) << k << "):";
        Editor e(instance, config.selector, config.forbidden, config.lower_bound);
        solved = e.edit(k,
                        [&](const std::vector<VertexPair> &edits) {
                            std::cout << Solution(instance, edits) << "\n";
                        },
                        [&](Cost pruned_k, Cost lb) { min_delta = std::min(min_delta, lb - pruned_k); deltas.push_back(lb - pruned_k); });
        std::sort(deltas.begin(), deltas.end(), [](Cost lhs, Cost rhs) { return lhs < rhs; });
        Cost delta10 = deltas[std::min(deltas.size() - 1, deltas.size() / 5)];
        std::cout << " min_delta = " << std::setw(10) << min_delta << ", deltas size = " << std::setw(10) << deltas.size() << ", delta10 = " << std::setw(10) << delta10 << "\n";
        if (!solved) k += delta10;
    } while (!solved);
}


int main() {

    const int multiplier = 100;
    const std::vector<std::string> paths {
        "../data/cost_matrix_component_nr_3_size_16_cutoff_10.0.metis",
        "../data/cost_matrix_component_nr_4_size_39_cutoff_10.0.metis",
        "../data/cost_matrix_component_nr_11_size_22_cutoff_10.0.metis",
        "./data/karate.graph"};
    auto instance = GraphIO::read_graph(paths[0], multiplier);

    auto selector = Options::Selector::LeastWeight;
    auto forbidden_type = Options::FSG::P4C4;
    auto lower_bound = Options::LB::No;


    std::vector<Solution> solutions;

    Editor editor(instance, selector, forbidden_type, lower_bound);
    auto solution_cb = [&](const std::vector<VertexPair> &edits) {
        Solution solution(instance, edits);
        std::cout << solution << "\n";
        solutions.push_back(solution);
    };
    auto pruning_cb = [](Cost k, Cost lb) { std::cout << "pruned: k=" << k << ", lb=" << lb << ", eps=" << lb - k << "\n"; };
    auto pruning_cb_2 = [](Cost, Cost) {};
    auto pruning_cb_3 = [](Cost k, Cost lb) { if (k != 0 || lb != 0) std::cout << "pruned: k=" << k << ", lb=" << lb << ", eps=" << lb - k << "\n"; };



    // 781 { {2, 6} {6, 8} {1, 4} {3, 8} {7, 8} {4, 11} {0, 3} }
    // 777 { {2, 6} {6, 8} {3, 8} {5, 8} {7, 8} {1, 8} {4, 11} {0, 3} }
    // 767 { {0, 3} {6, 8} {3, 8} {5, 8} {7, 8} {8, 9} {4, 10} }
    // 753 { {0, 3} {6, 8} {3, 8} {1, 8} {0, 7} {4, 11} }
    // 748 { {0, 3} {6, 8} {3, 8} {5, 8} {7, 8} {8, 9} {1, 8} {4, 11} }
    // 728 { {0, 2} {6, 8} {3, 8} {1, 8} {0, 7} {4, 11} }
    // 681 { {0, 3} {6, 8} {3, 8} {1, 4} {7, 8} {4, 11} }
    // 677 { {3, 6} {6, 8} {3, 8} {1, 8} {7, 8} {4, 11} }
    // 677 { {0, 3} {6, 8} {3, 8} {5, 8} {7, 8} {1, 8} {4, 11} }
    // 597 { {3, 6} {6, 8} {1, 8} {7, 8} {4, 11} }
    // 526 { {2, 6} {6, 8} {3, 8} {1, 8} {7, 8} {4, 11} {0, 3} }
    // 426 { {0, 3} {6, 8} {3, 8} {1, 8} {7, 8} {4, 11} }


    bool solved = editor.edit(6 * multiplier, solution_cb, pruning_cb);



    std::cout << (solved ? "instance solved" : "instance not solved") << "\n";

    std::sort(solutions.begin(), solutions.end(),
              [](const auto &lhs, const auto &rhs) { return lhs.cost > rhs.cost; });


    for (const auto& solution : solutions) {
        assert(is_solution_valid(instance.graph, solution, forbidden_type));
        std::cout << solution << "\n";
    }


    for (int seed = 0; seed < 1; ++seed) {
        FinderTests(seed).run();
        GraphTests(seed).run();
        SubgraphTests(seed).run();
        EditorTests(seed).run();
    }



    return 0;
}