//
// Created by jonas on 28.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FINDERTESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FINDERTESTS_H


#include "../graph/Graph.h"
#include "../Configuration.h"
#include "Tests.h"

#include "../interfaces/FinderI.h"
#include "../finder/NaiveC4P4.h"
#include "../finder/NaiveP3.h"
#include "../finder/CenterC4P4.h"


std::unique_ptr<FinderI> create_finder(const Graph& graph, Configuration::ForbiddenSubgraphs forbidden) {
    switch (forbidden) {
        case Configuration::ForbiddenSubgraphs::P3:
            return std::make_unique<Finder::NaiveP3>(graph);
        case Configuration::ForbiddenSubgraphs::P4C4:
            return std::make_unique<Finder::NaiveC4P4>(graph);
        default:
            abort();
    }
}


bool is_solution_valid(Graph &graph, const std::vector<VertexPair> &edits, Configuration::ForbiddenSubgraphs forbidden) {
    auto finder = create_finder(graph, forbidden);

    for (VertexPair uv : edits)
        graph.toggle_edge(uv);

    bool found_forbidden_subgraph = finder->find([&](const Subgraph &subgraph) { return true; });

    for (VertexPair uv : edits)
        graph.toggle_edge(uv);

    return !found_forbidden_subgraph;
}

std::vector<Subgraph> find_all_subgraphs(const Graph &graph, FinderI *finder) {
    std::vector<Subgraph> subgraphs;
    finder->find([&](Subgraph&& subgraph) { subgraphs.push_back(subgraph); return false; });
    return subgraphs;
}

std::vector<Subgraph> find_all_subgraphs(const Graph &graph, Configuration::ForbiddenSubgraphs forbidden) {
    auto finder = create_finder(graph, forbidden);
    return find_all_subgraphs(graph, finder.get());
}


class FinderTests {
    std::mt19937 gen;
public:
    explicit FinderTests(int seed=0) : gen(seed) {}

    void FinderFindsP4() {
        Graph G(4);
        G.set_edges({{0, 1}, {1, 2},{2, 3}});

        std::vector<Subgraph> expected{{0, 1, 2, 3}};
        auto actual = find_all_subgraphs(G, Configuration::ForbiddenSubgraphs::P4C4);
        expect("Finder recognizes P4", normalize(expected), normalize(actual));
    }

    void FinderFindsC4() {
        Graph G(4);
        G.set_edges({{0, 1}, {1, 2},{2, 3}, {3, 0}});

        std::vector<Subgraph> expected{{0, 1, 2, 3}};
        auto actual = find_all_subgraphs(G, Configuration::ForbiddenSubgraphs::P4C4);
        expect("Finder recognizes C4", normalize(expected), normalize(actual));
    }

    void FinderFindsP3() {
        Graph G(3);
        G.set_edges({{0, 1}, {1, 2}});

        std::vector<Subgraph> expected{{0, 1, 2}};
        auto actual = find_all_subgraphs(G, Configuration::ForbiddenSubgraphs::P3);
        expect("Finder recognizes P3", normalize(expected), normalize(actual));
    }

    void EditsSolveKarate()         {
        Graph G = GraphIO::read_graph("../data/karate.graph").graph;

        std::vector<VertexPair> edits {
            {0, 8}, {0, 16}, {0, 31}, /*{0, 33},*/
            /*{1, 3},*/ {1, 12}, {1, 30},
            {2, 8}, {2, 9}, {2, 27}, {2, 28}, {2, 32}, /*{2, 33},*/
            {4, 6}, {4, 10},
            {13, 33},
            {19, 33},
            {23, 25}, {13, 27},
            {24, 27},
            {26, 29},
            {28, 31},
            {31, 32}, {31, 33}
        };
        expect("Edits solve karate", true, is_solution_valid(G, edits, Configuration::ForbiddenSubgraphs::P4C4));
    }

    void NaiveP3_and_CenterC4P4_have_same_find_result() {
        Graph G = random_graph(10, 40, gen);

        auto F1 = std::make_unique<Finder::NaiveP3>(G);
        auto F2 = std::make_unique<Finder::CenterP3>(G);

        expect("Naive and Center P3 Finder are consistent", normalize(find_all_subgraphs(G, F1.get())), normalize(find_all_subgraphs(G, F2.get())));
    }

    void NaiveC4P4_and_CenterC4P4_have_same_find_result() {
        Graph G = random_graph(10, 40, gen);

        auto F1 = std::make_unique<Finder::NaiveC4P4>(G);
        auto F2 = std::make_unique<Finder::CenterC4P4>(G);

        expect("Naive and Center C4P4 Finder are consistent", normalize(find_all_subgraphs(G, F1.get())), normalize(find_all_subgraphs(G, F2.get())));
    }

    void run() {
        FinderFindsC4();
        FinderFindsP4();
        FinderFindsP3();
        EditsSolveKarate();
        NaiveP3_and_CenterC4P4_have_same_find_result();
        NaiveC4P4_and_CenterC4P4_have_same_find_result();
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FINDERTESTS_H
