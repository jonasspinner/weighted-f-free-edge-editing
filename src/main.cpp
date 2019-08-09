#include <iostream>
#include <random>

#include "graph/Graph.h"
#include "graph/GraphIO.h"
#include "Statistics.h"
#include "Editor.h"

#include "interfaces/FinderI.h"
#include "finder/NaiveP3.h"
#include "finder/CenterC4P4.h"

#include "lower_bound/ExplicitMWISSolver.h"
#include "lower_bound/NoLowerBound.h"

#include "selector/FirstEditable.h"
#include "selector/LeastWeight.h"
#include "tests/FinderTests.h"
#include "tests/GraphTests.h"
#include "tests/SubgraphTests.h"


class Solution {
public:
    Cost cost;
    std::vector<VertexPair> edits;

    // Solution(Cost cost, std::vector<VertexPair> edits) : cost(cost), edits(std::move(edits)) {}
    Solution(const VertexPairMap<Cost> &costs, std::vector<VertexPair> _edits) : cost(0), edits(std::move(_edits)) {
        for (VertexPair uv : edits) cost += costs[uv];
    }
    Solution(const Instance &instance, std::vector<VertexPair> edits) : Solution(instance.costs, std::move(edits)) {}

    friend std::ostream &operator<<(std::ostream &os, const Solution &solution) {
        os << solution.cost << " { ";
        for (auto uv : solution.edits) os << uv << " ";
        return os << "}";
    }
};

[[nodiscard]] bool is_solution_valid(const Graph &graph, const Solution &solution, Configuration::ForbiddenSubgraphs forbidden) {
    Graph graph_(graph);
    FinderI* finder;
    switch (forbidden) {
        case Configuration::ForbiddenSubgraphs::P3:
            finder = new Finder::NaiveP3(graph_);
            break;
        case Configuration::ForbiddenSubgraphs::P4C4:
            finder = new Finder::CenterC4P4(graph_);
            break;
    }

    for (VertexPair uv : solution.edits)
        graph_.toggle_edge(uv);

    bool found_forbidden_subgraph = finder->find([&](const Subgraph &subgraph) { return true; });

    return !found_forbidden_subgraph;
}

// s_b(s, n_b) * n_b >= s
// (s_b(s, n_b) + 1) * n_b < s


// (s - 1) // n_b

// 0                   s_b(1, 4) = 1
// 0 1                 s_b(2, 4) = 1
// 0 1 2               s_b(3, 4) = 1
// 0 1 2 3             s_b(4, 4) = 1
// 0 0 1 1 2           s_b(5, 4) = 2
// 0 0 1 1 2 2         s_b(6, 4) = 2
// 0 0 1 1 2 2 3       s_b(7, 4) = 2
// 0 0 1 1 2 2 3 3     s_b(8, 4) = 2
// 0 0 0 1 1 1 2 2 2   s_b(9, 4) = 3

// 0                   s_b(1, 3) = 1
// 0 1                 s_b(2, 3) = 1
// 0 1 2               s_b(3, 3) = 1
// 0 0 1 1             s_b(4, 3) = 2
// 0 0 1 1 2           s_b(5, 3) = 2
// 0 0 1 1 2 2         s_b(6, 3) = 2
// 0 0 0 1 1 1 2       s_b(7, 3) = 3
// 0 0 0 1 1 1 2 2     s_b(8, 3) = 3
// 0 0 0 1 1 1 2 2 2   s_b(9, 3) = 3

// 0                   s_b(1, 2) = 1
// 0 1                 s_b(1, 2) = 1
// 0 0 1               s_b(3, 2) = 2
// 0 0 1 1             s_b(4, 2) = 2
// 0 0 0 1 1           s_b(5, 2) = 3
// 0 0 0 1 1 1         s_b(6, 2) = 3
// 0 0 0 0 1 1 1       s_b(7, 2) = 4
// 0 0 0 0 1 1 1 1     s_b(8, 2) = 4
// 0 0 0 0 0 1 1 1 1   s_b(9, 2) = 5


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
                        [&](Cost k, Cost lb) { min_delta = std::min(min_delta, lb - k); deltas.push_back(lb -k); });
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

    auto selector = Configuration::SelectorOption::LeastWeight;
    auto forbidden = Configuration::ForbiddenSubgraphs::P4C4;
    auto lower_bound = Configuration::LowerBound::IteratedLocalSearch;


    std::vector<Solution> solutions;

    Editor editor(instance, selector, forbidden, lower_bound);
    auto solution_cb = [&](const std::vector<VertexPair> &edits) {
        Solution solution(instance, edits);
        std::cout << solution << "\n";
        solutions.push_back(solution);
    };
    auto pruning_cb = [](Cost k, Cost lb) { std::cout << "pruned: k=" << k << ", lb=" << lb << ", eps=" << lb - k << "\n"; };

    bool solved = editor.edit(60 * multiplier +  39 + 16, solution_cb, pruning_cb);

    std::cout << (solved ? "instance solved" : "instance not solved") << "\n";

    std::sort(solutions.begin(), solutions.end(),
              [](const auto &lhs, const auto &rhs) { return lhs.cost > rhs.cost; });


    for (const auto& solution : solutions) {
        assert(is_solution_valid(instance.graph, solution, forbidden));
        std::cout << solution << "\n";
    }


    for (int seed = 0; seed < 1; ++seed) {
        FinderTests(seed).run();
        GraphTests(seed).run();
        SubgraphTests(seed).run();
    }


    const Graph& G = instance.graph;

    std::cout << "vertices:";
    for (Vertex u : G.vertices())
        std::cout << " " << u;
    std::cout << "\n";

    std::cout << "vertex pairs:";
    for (VertexPair uv : G.vertexPairs())
        std::cout << " " << uv;
    std::cout << "\n";

    std::cout << "edges:";
    for (VertexPair uv : G.edges())
        std::cout << " " << uv;
    std::cout << "\n";

    std::cout << "neighbors(1):";
    for (Vertex u : G.neighbors(1))
        std::cout << " " << u;
    std::cout << "\n";

    /*
    auto c = std::make_unique<Finder::CenterP3>(G);
    auto n = std::make_unique<Finder::NaiveP3>(G);

    auto print_sg = [&G](auto sg) {
        for (Vertex u : sg) { std::cout << std::setw(2) << u << " "; }
        std::cout << G.has_edge({sg[0], sg[1]})
                  << G.has_edge({sg[0], sg[2]})
                  << G.has_edge({sg[1], sg[2]}) << "\n";
        return false;
    };

    auto get_count = [&G](FinderI *f, VertexPair uv) {
        int count = 0;
        f->find_near(uv, [&count](auto sg) { count++; return false; });
        return count;
    };

    std::vector<Subgraph> c_sg, n_sg;
    c->find([&](auto sg) { c_sg.push_back(sg); return false; });
    n->find([&](auto sg) { n_sg.push_back(sg); return false; });
    std::cout << c_sg.size() << " " << n_sg.size() << "\n";

    c_sg.clear(); n_sg.clear();
    c->find_near({0, 3}, [&](auto sg) { print_sg(sg); c_sg.push_back(sg); return false; });
    n->find_near({0, 3}, [&](auto sg) { print_sg(sg); n_sg.push_back(sg); return false; });
    std::cout << c_sg.size() << " " << n_sg.size() << "\n";

    std::mt19937 gen;
    std::uniform_int_distribution<Vertex> dist(0, G.size() - 2);
    for (int i = 0; i < 100; ++i) {
        Vertex u = dist(gen), v = dist(gen);
        v = (u != v) ? v : u + 1;
        assert(get_count(c.get(), {u, v}) == get_count(n.get(), {u, v}));
    }
    */
    // F_P3->find_near({0, 2}, print_sg);


    return 0;
}