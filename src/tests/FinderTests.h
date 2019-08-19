//
// Created by jonas on 28.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FINDERTESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FINDERTESTS_H


#include "../graph/Graph.h"
#include "../graph/GraphIO.h"
#include "../Configuration.h"
#include "../Solution.h"
#include "Tests.h"

#include "../interfaces/FinderI.h"
#include "../finder/NaiveC4P4.h"
#include "../finder/NaiveP3.h"
#include "../finder/CenterC4P4.h"
#include "../finder/CenterP3.h"
#include "../finder/Center.h"
#include "../finder/util.h"



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

std::vector<Subgraph> find_all_non_marked_subgraphs(FinderI &finder, const Graph& marked) {
    std::vector<Subgraph> subgraphs;
    finder.find(marked, [&](Subgraph&& subgraph) { subgraphs.push_back(subgraph); return false; });
    return subgraphs;
}


class FinderTests {
    std::mt19937 gen;
public:
    explicit FinderTests(int seed=0) : gen(seed) {}

    void EditsSolveKarate() {
        Instance instance = GraphIO::read_graph("./data/karate.graph");
        Graph G = instance.graph;

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
        expect("Edits solve karate", true, Solution(instance, edits).is_valid(instance, Options::FSG::P4C4));
    }

    template <typename A, typename B>
    void C4P4_Finders_are_consistent(const std::string& a_name, const std::string& b_name) {
        Graph G = random_graph(15, 200, gen);

        A a_finder(G); B b_finder(G);
        auto a_subgraphs = find_all_subgraphs(a_finder);
        auto b_subgraphs = find_all_subgraphs(b_finder);

        expect(a_name + " only produces C4P4", true, all_c4p4(G, a_subgraphs));
        expect(b_name + " only produces C4P4", true, all_c4p4(G, b_subgraphs));

        auto a_normalized = normalize(a_subgraphs);
        auto b_normalized = normalize(b_subgraphs);

        expect(a_name + " and " + b_name + " C4P4 Finder have same find output", a_normalized, b_normalized);

        a_normalized.erase(std::unique(a_normalized.begin(), a_normalized.end()), a_normalized.end());
        b_normalized.erase(std::unique(b_normalized.begin(), b_normalized.end()), b_normalized.end());
        expect(a_name + " and " + b_name + " C4P4 Finder have same find output ignoring duplicates", a_normalized, b_normalized);


        std::vector<std::vector<Subgraph>> a_near_subgraphs, b_near_subgraphs;
        for (VertexPair uv : G.vertexPairs()) {
            std::vector<Subgraph> a, b;

            a_finder.find_near(uv, [&](Subgraph &&subgraph){
                a.emplace_back(std::move(subgraph));
                return false;
            });
            a_near_subgraphs.emplace_back(normalize(a));

            b_finder.find_near(uv, [&](Subgraph &&subgraph){
                b.emplace_back(std::move(subgraph));
                return false;
            });
            b_near_subgraphs.emplace_back(normalize(b));
        }

        expect(a_name + " and " + b_name + " C4P4 Finder have same find_near output", a_near_subgraphs, b_near_subgraphs);

        for (auto &subgraphs : a_near_subgraphs) {
            subgraphs.erase(std::unique(subgraphs.begin(), subgraphs.end()), subgraphs.end());
        }
        for (auto &subgraphs : b_near_subgraphs) {
            subgraphs.erase(std::unique(subgraphs.begin(), subgraphs.end()), subgraphs.end());
        }

        expect(a_name + " and " + b_name + " C4P4 Finder have same find_near output ignoring duplicates", a_near_subgraphs, b_near_subgraphs);

    }

    template <typename Finder>
    void Finder_finds_C4(const std::string& name) {
        Graph C4(4);
        C4.set_edges({{0, 1}, {1, 2},{2, 3}, {3, 0}});

        {
            std::vector<Subgraph> expected{{0, 1, 2, 3}};

            Finder finder(C4);
            auto actual = find_all_subgraphs(finder);
            expect(name + " recognizes C4", normalize(expected), normalize(actual));
        }

        // 3-0-4
        // | |
        // 2-1-5
        Graph G(6);
        G.set_edges({{0, 1}, {1 , 2}, {2 , 3}, {3 , 0}, {0 , 4}, {1 , 5}});

        {
            std::vector<Subgraph> expected({{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 5}, {0, 1, 4, 5}, {0, 2, 3, 4}, {1, 2, 3, 5}});
            Finder finder(G);
            auto actual = find_all_subgraphs(finder);
            expect(name + " recognizes C4P4 in small graph", normalize(expected), normalize(actual));
        }

        {
            Graph marked(6);
            marked.set_edge({1, 2});
            std::vector<Subgraph> expected({{0, 1, 3, 5}, {0, 1, 4, 5}, {0, 2, 3, 4}});
            Finder finder(G);
            auto actual = find_all_non_marked_subgraphs(finder, marked);
            expect(name + " recognizes C4P4 in small graph with marked edges", normalize(expected), normalize(actual));
        }

        /*    4
              |
            3-0-5
             /|
          b-2-1-6
           /| |\
          a 9 8 7  */
        Graph G2(12);
        G2.set_edges({{0, 1}, {1, 2}, {2, 0}, {0, 3}, {0, 4}, {0, 5}, {1, 6}, {1, 7}, {1, 8}, {2, 9}, {2, 10}, {2, 11}});

        {
            Finder finder(G2);
            std::vector<Subgraph> expected({{3, 0, 1, 6}, {3, 0, 1, 7}, {3, 0, 1, 8}, {3, 0, 2, 9}, {3, 0, 2, 10}, {3, 0, 2, 11},
                                            {4, 0, 1, 6}, {4, 0, 1, 7}, {4, 0, 1, 8}, {4, 0, 2, 9}, {4, 0, 2, 10}, {4, 0, 2, 11},
                                            {5, 0, 1, 6}, {5, 0, 1, 7}, {5, 0, 1, 8}, {5, 0, 2, 9}, {5, 0, 2, 10}, {5, 0, 2, 11},
                                            {6, 1, 2, 9}, {6, 1, 2, 10}, {6, 1, 2, 11}, {7, 1, 2, 9}, {7, 1, 2, 10}, {7, 1, 2, 11}, {8, 1, 2, 9}, {8, 1, 2, 10}, {8, 1, 2, 11}});
            auto actual = find_all_subgraphs(finder);
            expect(name + " recognizes many P4 in small graph", normalize(expected), normalize(actual));
        }

        {
            Graph marked(12);
            marked.set_edge({0, 1});
            Finder finder(G2);
            std::vector<Subgraph> expected({{3, 0, 2, 9}, {3, 0, 2, 10}, {3, 0, 2, 11},
                                            {4, 0, 2, 9}, {4, 0, 2, 10}, {4, 0, 2, 11},
                                            {5, 0, 2, 9}, {5, 0, 2, 10}, {5, 0, 2, 11},
                                            {6, 1, 2, 9}, {6, 1, 2, 10}, {6, 1, 2, 11}, {7, 1, 2, 9}, {7, 1, 2, 10}, {7, 1, 2, 11}, {8, 1, 2, 9}, {8, 1, 2, 10}, {8, 1, 2, 11}});
            auto actual = find_all_non_marked_subgraphs(finder, marked);
            expect(name + " recognizes many P4 in small graph with marked edges", normalize(expected), normalize(actual));
        }
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


    template <typename A, typename B>
    void P3_Finders_are_consistent(const std::string& a_name, const std::string& b_name) {
        Graph G = random_graph(10, 40, gen);

        A a_finder(G); B b_finder(G);
        auto a_subgraphs = find_all_subgraphs(a_finder);
        auto b_subgraphs = find_all_subgraphs(b_finder);

        // expect(a_name + " only produces P3", true, all_p3(G, a_subgraphs));
        // expect(b_name + " only produces P3", true, all_p3(G, b_subgraphs));

        auto a_normalized = normalize(a_subgraphs);
        auto b_normalized = normalize(b_subgraphs);

        expect(a_name + " and " + b_name + " P3 Finder have same output", a_normalized, b_normalized);

        a_normalized.erase(std::unique(a_normalized.begin(), a_normalized.end()), a_normalized.end());
        b_normalized.erase(std::unique(b_normalized.begin(), b_normalized.end()), b_normalized.end());
        expect(a_name + " and " + b_name + " P3 Finder have same output ignoring duplicates", a_normalized, b_normalized);
    }

    template <typename Finder>
    void Finder_finds_P3(const std::string& name) {
        Graph G(3);
        G.set_edges({{0, 1}, {1, 2}});

        std::vector<Subgraph> expected{{0, 1, 2}};
        Finder finder(G);
        auto actual = find_all_subgraphs(finder);
        expect(name + " recognizes P3", normalize(expected), normalize(actual));
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

        Finder_finds_P3<Finder::CenterP3>("CenterP3");
        Finder_finds_P3<Finder::NaiveP3>("NaiveP3");

        P3_Finders_are_consistent<Finder::CenterP3, Finder::NaiveP3>("CenterP3", "NaiveP3");

        EditsSolveKarate();
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FINDERTESTS_H
