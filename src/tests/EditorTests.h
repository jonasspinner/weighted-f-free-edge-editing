//
// Created by jonas on 12.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_EDITORTESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_EDITORTESTS_H


#include "../Configuration.h"
#include "../graph/Graph.h"
#include "Tests.h"
#include "../Editor.h"

class EditorTests {

    std::mt19937 gen;
public:
    explicit EditorTests(int seed = 0) : gen(seed) {}

    void lower_bounds_have_same_output(const std::string &name_a, Options::LB lb_a, const std::string &name_b, Options::LB lb_b) {
        auto instance = GraphIO::read_graph("../data/cost_matrix_component_nr_3_size_16_cutoff_10.0.metis", 100);

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

     void run() {
         std::cout << "\nEditorTests"
                      "\n-----------" << std::endl;

         lower_bounds_have_same_output("No", Options::LB::No, "Greedy", Options::LB::Greedy);
         lower_bounds_have_same_output("No", Options::LB::No, "LocalSearch", Options::LB::LocalSearch);

     }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_EDITORTESTS_H
