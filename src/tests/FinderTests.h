//
// Created by jonas on 28.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FINDERTESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FINDERTESTS_H


#include "../graph/Graph.h"
#include "../graph/GraphIO.h"
#include "../Configuration.h"
#include "../Solution.h"
#include "test_utils.h"

#include "../interfaces/FinderI.h"
#include "../finder/NaiveC4P4.h"
#include "../finder/NaiveP3.h"
#include "../finder/CenterC4P4.h"
#include "../finder/CenterP3.h"
#include "../finder/Center.h"
#include "../finder/finder_utils.h"
#include "../Permutation.h"
#include "../finder/Endpoint.h"
#include "../finder/Naive.h"
#include "../finder/OuterP3.h"


bool is_p4(const Graph &graph, const Subgraph &subgraph) {
    const Subgraph& P = subgraph;
    assert(P.size() == 4);

    bool result = true;
    result &= graph.hasEdge({P[0], P[1]});
    result &= !graph.hasEdge({P[0], P[2]});
    result &= !graph.hasEdge({P[0], P[3]});
    result &= graph.hasEdge({P[1], P[2]});
    result &= !graph.hasEdge({P[1], P[3]});
    result &= graph.hasEdge({P[2], P[3]});
    return result;
}

bool is_c4(const Graph &graph, const Subgraph &subgraph) {
    const Subgraph& C = subgraph;
    assert(C.size() == 4);

    bool result = true;
    result &= graph.hasEdge({C[0], C[1]});
    result &= !graph.hasEdge({C[0], C[2]});
    result &= graph.hasEdge({C[0], C[3]});
    result &= graph.hasEdge({C[1], C[2]});
    result &= !graph.hasEdge({C[1], C[3]});
    result &= graph.hasEdge({C[2], C[3]});
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
    explicit FinderTests(int seed=0) : gen(static_cast<unsigned long>(seed)) {}

    static std::unique_ptr<FinderI> make_finder(const std::string &name, const Graph &graph) {
        if (name == "CenterRecC5P5") {
            return std::make_unique<Finder::CenterRecC5P5>(graph);
        } else if (name == "CenterRecC4P4") {
            return std::make_unique<Finder::CenterRecC4P4>(graph);
        } else if (name == "CenterRecP3") {
            return std::make_unique<Finder::CenterRecP3>(graph);
        } else if (name == "EndpointRecC5P5") {
            return std::make_unique<Finder::EndpointRecC5P5>(graph);
        } else if (name == "EndpointRecC4P4") {
            return std::make_unique<Finder::EndpointRecC4P4 >(graph);
        } else if (name == "EndpointRecP3") {
            return std::make_unique<Finder::EndpointRecP3>(graph);
        } else if (name == "CenterC4P4") {
            return std::make_unique<Finder::CenterC4P4>(graph);
        } else if (name == "CenterP3") {
            return std::make_unique<Finder::CenterP3>(graph);
        } else if (name == "NaiveC4P4") {
            return std::make_unique<Finder::NaiveC4P4>(graph);
        } else if (name == "NaiveP3") {
            return std::make_unique<Finder::NaiveP3>(graph);
        } else if (name == "NaiveRecC5P5") {
            return std::make_unique<Finder::NaiveRecC5P5>(graph);
        } else if (name == "NaiveRecC4P4") {
            return std::make_unique<Finder::NaiveRecC4P4>(graph);
        } else if (name == "NaiveRecP3") {
            return std::make_unique<Finder::NaiveRecP3>(graph);
        } else if (name == "OuterP3") {
            return std::make_unique<Finder::OuterP3>(graph);
        } else {
            std::cerr << "name = " << name << "\n";
            throw std::runtime_error("Finder name not valid.");
        }
    }

    void EditsSolveKarate() {
        try {
            auto instance = GraphIO::read_instance("../data/misc/karate.graph");
            Graph G = instance.graph;

            std::vector<VertexPair> edits {
                    {0, 8}, {0, 16}, {0, 31},
                    {1, 30},
                    {2, 8}, {2, 9}, {2, 27}, {2, 28}, {2, 32},
                    {3, 12},
                    {4, 5}, {4, 10},
                    {13, 33},
                    {19, 33},
                    {23, 25}, {23, 27},
                    {26, 29},
                    {27, 33},
                    {28, 31},
                    {31, 32}, {31, 33}
            };
            expect("Edits solve karate", true, Solution(instance, edits).is_valid(instance, Options::FSG::C4P4));

        } catch (const std::runtime_error &e) {
            std::cerr << e.what() << "\n";
        }
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
        C4.setEdges({{0, 1},
                     {1, 2},
                     {2, 3},
                     {3, 0}});

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
        G.setEdges({{0, 1},
                    {1, 2},
                    {2, 3},
                    {3, 0},
                    {0, 4},
                    {1, 5}});

        {
            std::vector<Subgraph> expected({{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 5}, {0, 1, 4, 5}, {0, 2, 3, 4}, {1, 2, 3, 5}});
            Finder finder(G);
            auto actual = find_all_subgraphs(finder);
            expect(name + " recognizes C4P4 in small graph", normalize(expected), normalize(actual));
        }

        {
            Graph marked(6);
            marked.setEdge({1, 2});
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
        G2.setEdges({{0, 1},
                     {1, 2},
                     {2, 0},
                     {0, 3},
                     {0, 4},
                     {0, 5},
                     {1, 6},
                     {1, 7},
                     {1, 8},
                     {2, 9},
                     {2, 10},
                     {2, 11}});

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
            marked.setEdge({0, 1});
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
        G.setEdges({{0, 1},
                    {1, 2},
                    {2, 3}});

        std::vector<Subgraph> expected{{0, 1, 2, 3}};
        Finder finder(G);
        auto actual = find_all_subgraphs(finder);
        expect(name + " recognizes P4", normalize(expected), normalize(actual));
    }

    template <typename Finder>
    void Finder_finds_P3(const std::string& name) {
        Graph G(3);
        G.setEdges({{0, 1},
                    {1, 2}});

        std::vector<Subgraph> expected{{0, 1, 2}};
        Finder finder(G);
        auto actual = find_all_subgraphs(finder);
        expect(name + " recognizes P3", normalize(expected), normalize(actual));
    }


    template <class Finder>
    void Finder_is_seed_independent(const std::string &name, const std::vector<int> &seeds) {
        Graph G = random_graph(10, 40, gen);

        std::vector<std::pair<int, std::vector<Subgraph>>> results;
        for (auto seed : seeds) {
            auto P = Permutation(G.size(), seed);
            auto P_r = P.reverse();

            auto G_p = P[G];

            std::vector<Subgraph> subgraphs;
            Finder finder(G_p);
            finder.find([&](Subgraph &&subgraph) {
                subgraphs.push_back(P_r[subgraph]);
                return false;
            });

            results.emplace_back(seed, std::move(subgraphs));
        }

        for (size_t i = 0; i < seeds.size(); ++i) {
            for (size_t j = i + 1; j < seeds.size(); ++j) {
                const auto &[seed_i, subgraphs_i] = results[i];
                const auto &[seed_j, subgraphs_j] = results[j];
                std::stringstream ss;
                ss << "Finder " << name << " finds the same subgraphs with seed " << seed_i << " and " << seed_j;
                expect(ss.str(), normalize(subgraphs_i), normalize(subgraphs_j));
            }
        }
    }

    void finders_have_same_output(const std::vector<std::string> &finders, const std::vector<int> &seeds = {0}) {
        Graph G_original = random_graph(20, 200, gen);
        Graph F_original = random_graph(20, 50, gen);

        std::vector<std::tuple<std::string, int,
                    std::vector<Subgraph>, std::vector<Subgraph>,
                    std::vector<std::vector<Subgraph>>, std::vector<std::vector<Subgraph>>>> results;

        for (auto seed : seeds) {
            Permutation P(G_original.size(), seed);
            auto G = P[G_original];
            auto F = P[F_original];
            for (const auto &name : finders) {
                auto finder = make_finder(name, G);

                std::vector<Subgraph> find;
                finder->find([&](Subgraph &&subgraph) {
                    find.push_back(std::move(subgraph));
                    return false;
                });

                std::vector<Subgraph> find_forbidden;
                finder->find(F, [&](Subgraph &&subgraph) {
                    find_forbidden.push_back(std::move(subgraph));
                    return false;
                });

                std::vector<std::vector<Subgraph>> find_near;
                for (VertexPair uv : G.vertexPairs()) {
                    std::vector<Subgraph> find_near_uv;
                    finder->find_near(uv, [&](Subgraph &&subgraph) {
                        find_near_uv.push_back(std::move(subgraph));
                        return false;
                    });
                    find_near.push_back(normalize(find_near_uv));
                }

                std::vector<std::vector<Subgraph>> find_near_forbidden;
                for (VertexPair uv : G.vertexPairs()) {
                    std::vector<Subgraph> find_near_uv;
                    finder->find_near(uv, F, [&](Subgraph &&subgraph) {
                        find_near_uv.push_back(std::move(subgraph));
                        return false;
                    });
                    find_near_forbidden.push_back(normalize(find_near_uv));
                }

                results.emplace_back(name, seed, normalize(find), normalize(find_forbidden), find_near, find_near_forbidden);
            }
        }

        for (size_t i = 0; i < results.size(); ++i) {
            for (size_t j = i + 1; j < results.size(); ++j) {
                const auto &[name_i, seed_i, find_i, find_forbidden_i, find_near_i, find_near_forbidden_i] = results[i];
                const auto &[name_j, seed_j, find_j, find_forbidden_j, find_near_j, find_near_forbidden_j] = results[j];
                if (seed_i != seed_j) continue;
                {
                    std::stringstream ss;
                    ss << name_i << " and " << name_j << "have the same find output, seed=" << seed_i;
                    expect(ss.str(), find_i, find_j);
                }
                {
                    std::stringstream ss;
                    ss << name_i << " and " << name_j << "have the same find output with forbidden vertex pairs, seed=" << seed_i;
                    expect(ss.str(), find_forbidden_i, find_forbidden_j);
                }
                {
                    std::stringstream ss;
                    ss << name_i << " and " << name_j << "have the same find_near output, seed=" << seed_i;
                    expect(ss.str(), find_near_i, find_near_j);
                }
                {
                    std::stringstream ss;
                    ss << name_i << " and " << name_j << "have the same find_near output with forbidden vertex pairs, seed=" << seed_i;
                    expect(ss.str(), find_near_forbidden_i, find_near_forbidden_j);
                }
            }
        }
    }

    void run() {
        std::cout << "\nFinderTests"
                     "\n-----------" << std::endl;

        Finder_finds_P3<Finder::NaiveP3>("NaiveP3");
        Finder_finds_P3<Finder::CenterP3>("CenterP3");
        Finder_finds_P3<Finder::CenterRecP3>("CenterRecP3");
        Finder_finds_P3<Finder::EndpointRecP3>("EndpointRecP3");
        Finder_finds_P3<Finder::NaiveRecP3>("NaiveRecP3");
        Finder_finds_P3<Finder::NaiveRecP3>("OuterP3");
        finders_have_same_output({"NaiveP3", "CenterP3", "CenterRecP3", "EndpointRecP3", "NaiveRecP3", "OuterP3"}, {0, 1});

        Finder_finds_C4<Finder::NaiveC4P4>("NaiveC4P4");
        Finder_finds_C4<Finder::CenterC4P4>("CenterC4P4");
        Finder_finds_C4<Finder::CenterRecC4P4>("CenterRecC4P4");
        Finder_finds_C4<Finder::EndpointRecC4P4>("EndpointRecC4P4");
        Finder_finds_C4<Finder::NaiveRecC4P4>("NaiveRecC4P4");
        Finder_finds_P4<Finder::NaiveC4P4>("NaiveC4P4");
        Finder_finds_P4<Finder::CenterC4P4>("CenterC4P4");
        Finder_finds_P4<Finder::CenterRecC4P4>("CenterRecC4P4");
        Finder_finds_P4<Finder::EndpointRecC4P4>("EndpointRecC4P4");
        Finder_finds_P4<Finder::NaiveRecC4P4>("NaiveRecC4P4");
        finders_have_same_output({"NaiveC4P4", "CenterC4P4", "CenterRecC4P4", "EndpointRecC4P4", "NaiveRecC4P4"}, {0, 1});

        finders_have_same_output({"CenterRecC5P5", "EndpointRecC5P5", "NaiveRecC5P5"}, {0, 1});

        EditsSolveKarate();
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FINDERTESTS_H
