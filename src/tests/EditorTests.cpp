//
// Created by jonas on 12.08.19.
//


#include "EditorTests.h"
#include "../graph/GraphIO.h"
#include "../Solution.h"
#include "../editor/Editor.h"
#include "test_utils.h"


/**
 * Tests if all parameter combinations result in the same output.
 * The output is filtered by only considering inclusion minimal solutions.
 *
 * @param fsg
 * @param selectors
 * @param lower_bounds
 * @param seeds
 */
void EditorTests::configurations_have_same_output(Options::FSG fsg, const std::vector<Options::Selector> &selectors,
                                     const std::vector<Options::LB> &lower_bounds, const std::vector<int> &seeds,
                                     double multiplier) {
    std::cout << "\n#\tconfigurations_have_same_output\n";

    auto orig_instance = GraphIO::read_instance(instance_path, multiplier);
    std::vector<std::tuple<Options::Selector, Options::LB, int, std::vector<Solution>>> results;

    for (auto seed : seeds) {

        Permutation P(orig_instance.graph.size(), seed);
        Permutation P_r = P.reverse();

        auto instance = P[orig_instance];


        for (auto selector : selectors) {
            for (auto lb : lower_bounds) {
                auto config = Configuration(fsg, multiplier, Options::SolverType::FPT, selector, lb);
                config.input_path = instance_path;
                config.seed = seed;
                std::vector<Solution> solutions;

                Editor editor(instance.copy(), config);

                editor.edit(12 * multiplier, [&](const std::vector<VertexPair> &edits) {
                    solutions.emplace_back(orig_instance, P_r[edits]);
                }, [](Cost, Cost) {});

                results.emplace_back(selector, lb, seed, std::move(solutions));
            }
        }
    }


    for (auto &[_0, _1, _2, solutions] : results)
        Solution::filter_inclusion_minimal(solutions);


    // Compare pairwise
    for (size_t i = 0; i < results.size(); ++i) {
        for (size_t j = i + 1; j < results.size(); ++j) {
            const auto &[selector_i, lb_i, seed_i, solutions_i] = results[i];
            const auto &[selector_j, lb_j, seed_j, solutions_j] = results[j];

            std::stringstream name;
            name << "Editor(" << selector_i << ", " << fsg << ", " << lb_i << ") seed=" << seed_i << " and "
                 << "Editor(" << selector_j << ", " << fsg << ", " << lb_j << ") seed=" << seed_j << " "
                 << "have the same solutions";
            expect(name.str(), solutions_i, solutions_j);
        }
    }
}


void EditorTests::same_output_for_small_zero_cost_instance(const std::vector<Options::Selector> &selectors,
                                                           const std::vector<Options::LB> &lower_bounds,
                                                           const std::vector<int> &seeds) {
    std::cout << "\n#\tsame_output_for_small_zero_cost_instance\n";
    std::vector<std::tuple<Options::Selector, Options::LB, int, std::vector<Solution>>> results;

    auto fsg = Options::FSG::C4P4;
    auto graph = Graph::make_path_graph(5);
    VertexPairMap<Cost> costs(graph.size());
    for (auto uv : graph.vertexPairs()) {
        costs[uv] = 0;
    }
    Instance orig_instance{"test-instance", std::move(graph), costs};

    for (auto seed : seeds) {

        Permutation P(orig_instance.graph.size(), seed);
        Permutation P_r = P.reverse();

        auto instance = P[orig_instance];

        for (auto selector : selectors) {
            for (auto lb : lower_bounds) {
                auto config = Configuration(fsg, 1, Options::SolverType::FPT, selector, lb);
                std::vector<Solution> solutions;

                Editor editor(instance.copy(), config);

                editor.edit(1, [&](const std::vector<VertexPair> &edits) {
                    solutions.emplace_back(orig_instance, P_r[edits]);
                }, [](Cost, Cost) {});

                results.emplace_back(selector, lb, seed, std::move(solutions));
            }
        }
    }

    for (auto &[_0, _1, _2, solutions] : results)
        Solution::filter_inclusion_minimal(solutions);

    // Compare pairwise
    for (size_t i = 0; i < results.size(); ++i) {
        for (size_t j = i + 1; j < results.size(); ++j) {
            const auto &[selector_i, lb_i, seed_i, solutions_i] = results[i];
            const auto &[selector_j, lb_j, seed_j, solutions_j] = results[j];

            std::stringstream name;
            name << "Editor(" << selector_i << ", " << fsg << ", " << lb_i << ") seed=" << seed_i << " and "
                 << "Editor(" << selector_j << ", " << fsg << ", " << lb_j << ") seed=" << seed_j << " "
                 << "have the same solutions";
            expect(name.str(), solutions_i, solutions_j);
        }
    }
}


/**
 * Tests if the Editor class has output independent of the chosen seed.
 * The output is filtered by only considering inclusion minimal solutions.
 *
 * @param seeds
 */
void EditorTests::output_is_independent_of_seed(const std::vector<int> &seeds,
        Options::Selector selector = Options::Selector::FirstFound,
        Options::LB lb = Options::LB::SortedGreedy) {
    std::cout << "\n#\toutput_is_independent_of_seed\n";

    auto orig_instance = GraphIO::read_instance(instance_path, 100);

    std::vector<std::tuple<int, std::vector<Solution>>> results;

    for (int seed : seeds) {
        Permutation P(orig_instance.graph.size(), seed);
        Permutation P_r = P.reverse();

        auto instance = P[orig_instance];

        std::vector<Solution> solutions;

        auto config = Configuration(Options::FSG::C4P4, -1, Options::SolverType::FPT, selector, lb);

        Editor editor(instance.copy(), config);
        editor.edit(1100, [&](const std::vector<VertexPair> &edits) {
            solutions.emplace_back(orig_instance, P_r[edits]);
        }, [](Cost, Cost) {});

        results.emplace_back(seed, std::move(solutions));
    }


    for (auto &[_, solutions] : results)
        Solution::filter_inclusion_minimal(solutions);


    // Compare pairwise
    for (size_t i = 0; i < seeds.size(); ++i) {
        for (size_t j = i + 1; j < seeds.size(); ++j) {
            const auto &[seed_i, solutions_i] = results[i];
            const auto &[seed_j, solutions_j] = results[j];
            std::stringstream name;
            name << "seed = " << seed_i << " and "
                 << "seed = " << seed_j << " "
                 << "have the same solutions";
            expect(name.str(), solutions_i, solutions_j);
        }
    }
}

void EditorTests::run() {
    std::cout << "\nEditorTests"
                 "\n-----------" << std::endl;

    using Options::Selector;
    using Options::LB;
    using Options::FSG;

    // auto selectors = {Selector::FirstFound, Selector::LeastWeight, Selector::MostMarkedPairs, Selector::MostAdjacentSubgraphs};
    // auto lower_bounds = {LB::LocalSearch, LB::Trivial, LB::Greedy, LB::SortedGreedy, LB::LPRelaxation, LB::NPS_MWIS_Solver};

    auto selectors = {Selector::MostAdjacentSubgraphs};
    auto lower_bounds = {LB::WeightedPackingLocalSearch, LB::SortedGreedy, LB::GreedyWeightedPacking};

    configurations_have_same_output(FSG::C4P4, selectors, lower_bounds, {0, 1}, 100);
    configurations_have_same_output(FSG::C4P4, selectors, lower_bounds, {0, 1}, 1);

    //TODO: test case P4 with 0 cost vertex pairs.
    same_output_for_small_zero_cost_instance(selectors, lower_bounds, {0, 1});
    //configurations_have_same_output(FSG::C4P4, selectors, lower_bounds, {0}, 0);

    output_is_independent_of_seed({0, 1, 2, 3});
}

