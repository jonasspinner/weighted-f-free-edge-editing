//
// Created by jonas on 12.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_EDITORTESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_EDITORTESTS_H


#include "../Configuration.h"
#include "../graph/Graph.h"
#include "Tests.h"
#include "../Editor.h"
#include "../Permutation.h"

class EditorTests {

    std::mt19937 gen;
    std::string instance_path = "../data/cost_matrix_component_nr_3_size_16_cutoff_10.0.metis";
public:
    explicit EditorTests(int seed = 0) : gen(seed) {}

    void lower_bounds_have_same_output(const std::string &name_a, Options::LB lb_a, const std::string &name_b, Options::LB lb_b) {
        auto instance = GraphIO::read_graph(instance_path, 100);

        auto selector = Options::Selector::FirstEditable;
        auto forbidden = Options::FSG::P4C4;


        std::vector<std::vector<VertexPair>> edits_a, edits_b;
        {
            Editor editor_a(instance, selector, forbidden, lb_a);
            editor_a.edit(550, [&](const std::vector<VertexPair>& edits) { edits_a.push_back(edits); }, [](Cost, Cost) {});
        }

        {
            Editor editor_b(instance, selector, forbidden, lb_b);
            editor_b.edit(550, [&](const std::vector<VertexPair>& edits) { edits_b.push_back(edits); }, [](Cost, Cost) {});
        }

        expect(name_a + " and " + name_b + " have the same solutions", normalize(edits_a), normalize(edits_b));
    }

    void configurations_have_same_output(Options::FSG fsg, const std::vector<Options::Selector>& selectors, const std::vector<Options::LB>& lower_bounds, int seed = 0) {

        auto orig_instance = GraphIO::read_graph(instance_path, 100);

        Permutation P(orig_instance.graph.size(), seed);
        Permutation P_r = P.reverse();

        auto instance = P[orig_instance];

        std::vector<std::tuple<Options::Selector, Options::FSG, Options::LB, std::vector<Solution>>> results;

        for (auto selector : selectors) {
            for (auto lb : lower_bounds) {
                std::vector<Solution> solutions;

                Editor editor(instance, selector, fsg, lb);
                editor.edit(600, [&](const std::vector<VertexPair> &edits){
                    Cost cost = 0;
                    for (VertexPair uv : edits)
                        cost += instance.costs[uv];

                    solutions.emplace_back(instance, edits);
                }, [](Cost, Cost){});

                std::sort(solutions.begin(), solutions.end());
                results.emplace_back(selector, fsg, lb, std::move(solutions));
            }
        }

        for (size_t i = 0; i < results.size(); ++i) {
            for (size_t j = i + 1; j < results.size(); ++j) {
                const auto &[selector_i, fsg_i, lb_i, solutions_i] = results[i];
                const auto &[selector_j, fsg_j, lb_j, solutions_j] = results[i];
                expect("a b", solutions_i, solutions_j);
            }
        }
    }

     void run() {
         std::cout << "\nEditorTests"
                      "\n-----------" << std::endl;

         lower_bounds_have_same_output("No", Options::LB::No, "Greedy", Options::LB::Greedy);
         lower_bounds_have_same_output("No", Options::LB::No, "LocalSearch", Options::LB::LocalSearch);

         configurations_have_same_output(Options::FSG::P4C4, {Options::Selector::FirstEditable}, {Options::LB::Greedy});
     }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_EDITORTESTS_H
