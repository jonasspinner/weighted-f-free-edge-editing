//
// Created by jonas on 12.08.19.
//


#include "EditorTests.h"
#include "../graph/GraphIO.h"
#include "../Solution.h"
#include "../Permutation.h"
#include "../Editor.h"
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

    auto orig_instance = GraphIO::read_graph(instance_path, multiplier);
    std::vector<std::tuple<Options::Selector, Options::FSG, Options::LB, int, std::vector<Solution>>> results;

    for (auto seed : seeds) {

        Permutation P(orig_instance.graph.size(), seed);
        Permutation P_r = P.reverse();

        auto instance = P[orig_instance];


        for (auto selector : selectors) {
            for (auto lb : lower_bounds) {
                std::vector<Solution> solutions;

                Editor editor(instance, selector, fsg, lb);
                editor.initialize(10 * multiplier);
                editor.edit(10 * multiplier, [&](const std::vector<VertexPair> &edits) {
                    solutions.emplace_back(orig_instance, P_r[edits]);
                }, [](Cost, Cost) {});

                results.emplace_back(selector, fsg, lb, seed, std::move(solutions));
            }
        }
    }


    for (auto &[_0, _1, _2, _3, solutions] : results)
        Solution::filter_inclusion_minimal(solutions);


    // Compare pairwise
    for (size_t i = 0; i < results.size(); ++i) {
        for (size_t j = i + 1; j < results.size(); ++j) {
            const auto &[selector_i, fsg_i, lb_i, seed_i, solutions_i] = results[i];
            const auto &[selector_j, fsg_j, lb_j, seed_j, solutions_j] = results[j];

            std::stringstream name;
            name << "Editor(" << selector_i << ", " << fsg_i << ", " << lb_i << ") seed=" << seed_i << " and "
                 << "Editor(" << selector_j << ", " << fsg_j << ", " << lb_j << ") seed=" << seed_j << " "
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
void EditorTests::output_is_independent_of_seed(const std::vector<int> &seeds) {
    auto orig_instance = GraphIO::read_graph(instance_path, 100);

    std::vector<std::tuple<int, std::vector<Solution>>> results;

    for (int seed : seeds) {
        Permutation P(orig_instance.graph.size(), seed);
        Permutation P_r = P.reverse();

        auto instance = P[orig_instance];

        std::vector<Solution> solutions;

        Editor editor(instance, Options::Selector::FirstEditable, Options::FSG::P4C4, Options::LB::Greedy);
        editor.initialize(1100);
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

    auto all_selectors = {Selector::FirstEditable, Selector::LeastWeight, Selector::MostMarked};
    auto all_lower_bounds = {LB::No, LB::Greedy, LB::LocalSearch, LB::LinearProgram};

    configurations_have_same_output(FSG::P4C4, all_selectors, all_lower_bounds, {0, 1}, 100);
    configurations_have_same_output(FSG::P4C4, all_selectors, all_lower_bounds, {0, 1}, 1);

    output_is_independent_of_seed({0, 1, 2, 3, 4});
}