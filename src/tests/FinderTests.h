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
#include "../finder/Center.h"


std::unique_ptr<FinderI> create_finder(const Graph& graph, Configuration::ForbiddenSubgraphs forbidden) {
    switch (forbidden) {
        case Configuration::ForbiddenSubgraphs::P3:
            return std::make_unique<Finder::CenterP3>(graph);
        case Configuration::ForbiddenSubgraphs::P4C4:
            return std::make_unique<Finder::CenterC4P4>(graph);
        default:
            abort();
    }
}


bool is_solution_valid(Graph &graph, const std::vector<VertexPair> &edits, Configuration::ForbiddenSubgraphs forbidden) {
    auto finder = create_finder(graph, forbidden);

    for (VertexPair uv : edits)
        graph.toggle_edge(uv);

    bool found_forbidden_subgraph = finder->find([&](const Subgraph & /* subgraph */) { return true; });

    for (VertexPair uv : edits)
        graph.toggle_edge(uv);

    return !found_forbidden_subgraph;
}

bool is_p4(const Graph &graph, const Subgraph &subgraph) {
    const Subgraph& P = subgraph;
    assert(P.size() == 4);

    bool result = true;
    result &=  graph.has_edge({P[0], P[1]});
    result &= !graph.has_edge({P[0], P[2]});
    result &= !graph.has_edge({P[0], P[3]});
    result &=  graph.has_edge({P[1], P[2]});
    result &= !graph.has_edge({P[1], P[3]});
    result &=  graph.has_edge({P[2], P[3]});
    return result;
}

bool is_c4(const Graph &graph, const Subgraph &subgraph) {
    const Subgraph& C = subgraph;
    assert(C.size() == 4);

    bool result = true;
    result &=  graph.has_edge({C[0], C[1]});
    result &= !graph.has_edge({C[0], C[2]});
    result &=  graph.has_edge({C[0], C[3]});
    result &=  graph.has_edge({C[1], C[2]});
    result &= !graph.has_edge({C[1], C[3]});
    result &=  graph.has_edge({C[2], C[3]});
    return result;
}

bool all_c4p4(const Graph &graph, const std::vector<Subgraph> &subgraphs) {
    bool result = true;
    for (const auto& subgraph : subgraphs) {
        result &= is_p4(graph, subgraph) || is_c4(graph, subgraph);
    }
    return result;
}


std::vector<Subgraph> find_all_subgraphs(FinderI &finder) {
    std::vector<Subgraph> subgraphs;
    finder.find([&](Subgraph&& subgraph) { subgraphs.push_back(subgraph); return false; });
    return subgraphs;
}

std::vector<Subgraph> find_all_subgraphs(const Graph &graph, Configuration::ForbiddenSubgraphs forbidden) {
    auto finder = create_finder(graph, forbidden);
    return find_all_subgraphs(*finder);
}


class FinderTests {
    std::mt19937 gen;
public:
    explicit FinderTests(int seed=0) : gen(seed) {}

    void FinderFindsP3() {
        Graph G(3);
        G.set_edges({{0, 1}, {1, 2}});

        std::vector<Subgraph> expected{{0, 1, 2}};
        Finder::CenterP3 finder(G);
        auto actual = find_all_subgraphs(finder);
        expect("CenterP3 recognizes P3", normalize(expected), normalize(actual));
    }

    void EditsSolveKarate()         {
        Graph G = GraphIO::read_graph("./data/karate.graph").graph;

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

    template <typename A, typename B>
    void C4P4_Finders_are_consistent(const std::string& a_name, const std::string& b_name) {
        Graph G = random_graph(10, 40, gen);

        A a_finder(G); B b_finder(G);
        auto a_subgraphs = find_all_subgraphs(a_finder);
        auto b_subgraphs = find_all_subgraphs(b_finder);

        expect(a_name + " only produces C4P4", true, all_c4p4(G, a_subgraphs));
        expect(b_name + " only produces C4P4", true, all_c4p4(G, b_subgraphs));

        auto a_normalized = normalize(a_subgraphs);
        auto b_normalized = normalize(b_subgraphs);

        expect(a_name + " and " + b_name + " C4P4 Finder have same output", a_normalized, b_normalized);

        a_normalized.erase(std::unique(a_normalized.begin(), a_normalized.end()), a_normalized.end());
        b_normalized.erase(std::unique(b_normalized.begin(), b_normalized.end()), b_normalized.end());
        expect(a_name + " and " + b_name + " C4P4 Finder have same output ignoring duplicates", a_normalized, b_normalized);
    }

    template <typename Finder>
    void Finder_finds_C4(const std::string& name) {
        Graph G(4);
        G.set_edges({{0, 1}, {1, 2},{2, 3}, {3, 0}});

        std::vector<Subgraph> expected{{0, 1, 2, 3}};
        Finder finder(G);
        auto actual = find_all_subgraphs(finder);
        expect(name + " recognizes C4", normalize(expected), normalize(actual));
    }

    template <typename Finder>
    void Finder_finds_P4(const std::string& name) {
        Graph G(4);
        G.set_edges({{0, 1}, {1, 2},{2, 3}});

        std::vector<Subgraph> expected{{0, 1, 2, 3}};
        Finder finder(G);
        auto actual = find_all_subgraphs(finder);
        expect(name + " recognizes P4", normalize(expected), normalize(actual));
    }

    void run() {
        std::cout << "\nFinderTests"
                     "\n-----------" << std::endl;

        Finder_finds_C4<Finder::NaiveC4P4>("NaiveC4P4");
        Finder_finds_C4<Finder::CenterC4P4>("CenterC4P4");
        Finder_finds_C4<CenterRecC4P4>("CenterRecC4P4");

        Finder_finds_P4<Finder::NaiveC4P4>("NaiveC4P4");
        Finder_finds_P4<Finder::CenterC4P4>("CenterC4P4");
        Finder_finds_P4<CenterRecC4P4>("CenterRecC4P4");

        C4P4_Finders_are_consistent<Finder::NaiveC4P4, Finder::CenterC4P4>("NaiveC4P4", "CenterC4P4");
        C4P4_Finders_are_consistent<Finder::NaiveC4P4, CenterRecC4P4>("NaiveC4P4", "CenterRecC4P4");
        C4P4_Finders_are_consistent<Finder::CenterC4P4, CenterRecC4P4>("CenterC4P4", "CenterRecC4P4");

        FinderFindsP3();
        EditsSolveKarate();
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FINDERTESTS_H
