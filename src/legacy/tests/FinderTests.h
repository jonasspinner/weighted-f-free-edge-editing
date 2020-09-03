//
// Created by jonas on 28.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FINDERTESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FINDERTESTS_H


#include "../graph/Graph.h"
#include "../graph/GraphIO.h"
#include "../Configuration.h"
#include "../Solution.h"m
#include "test_utils.h"

#include "../finder/FinderI.h"
#include "../finder/NaiveC4P4.h"
#include "../finder/NaiveP3.h"
#include "../finder/CenterC4P4.h"
#include "../finder/CenterP3.h"
#include "../finder/Center.h"
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


std::vector<Subgraph> find_all_subgraphs(FinderI &finder, const Graph& graph) {
    std::vector<Subgraph> subgraphs;
    finder.find(graph, [&](Subgraph&& subgraph) { subgraphs.push_back(subgraph); return false; });
    return subgraphs;
}

std::vector<Subgraph> find_all_subgraphs_with_duplicates(FinderI &finder, const Graph& graph, const Graph& forbidden) {
    std::vector<Subgraph> subgraphs;
    finder.find_with_duplicates(graph, forbidden,
            [&](Subgraph&& subgraph) { subgraphs.push_back(subgraph); return false; });
    return subgraphs;
}

std::vector<Subgraph> find_all_non_forbidden_subgraphs(FinderI &finder, const Graph& graph, const Graph& forbidden) {
    std::vector<Subgraph> subgraphs;
    finder.find(graph, forbidden, [&](Subgraph&& subgraph) { subgraphs.push_back(subgraph); return false; });
    return subgraphs;
}


class FinderTests {
    std::mt19937 gen;
public:
    explicit FinderTests(int seed=0) : gen(static_cast<unsigned long>(seed)) {}

    static std::unique_ptr<FinderI> make_finder(const std::string &name) {
        if (name == "CenterRecC5P5") {
            return std::make_unique<Finder::CenterRecC5P5>();
        } else if (name == "CenterRecC4P4") {
            return std::make_unique<Finder::CenterRecC4P4>();
        } else if (name == "CenterRecP3") {
            return std::make_unique<Finder::CenterRecP3>();
        } else if (name == "EndpointRecC5P5") {
            return std::make_unique<Finder::EndpointRecC5P5>();
        } else if (name == "EndpointRecC4P4") {
            return std::make_unique<Finder::EndpointRecC4P4 >();
        } else if (name == "EndpointRecP3") {
            return std::make_unique<Finder::EndpointRecP3>();
        } else if (name == "CenterC4P4") {
            return std::make_unique<Finder::CenterC4P4>();
        } else if (name == "CenterP3") {
            return std::make_unique<Finder::CenterP3>();
        } else if (name == "NaiveC4P4") {
            return std::make_unique<Finder::NaiveC4P4>();
        } else if (name == "NaiveP3") {
            return std::make_unique<Finder::NaiveP3>();
        } else if (name == "NaiveRecC5P5") {
            return std::make_unique<Finder::NaiveRecC5P5>();
        } else if (name == "NaiveRecC4P4") {
            return std::make_unique<Finder::NaiveRecC4P4>();
        } else if (name == "NaiveRecP3") {
            return std::make_unique<Finder::NaiveRecP3>();
        } else if (name == "OuterP3") {
            return std::make_unique<Finder::OuterP3>();
        } else {
            std::cerr << "name = " << name << "\n";
            throw std::runtime_error("Finder name not valid.");
        }
    }

    void FinderC4P4_karate_graph_is_solved_with_21_edits() {
        try {
            auto instance = GraphIO::read_instance("../data/misc/karate.graph");
            Graph G = instance.graph.copy();

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
        auto C4 = Graph::make_cycle_graph(4);

        {
            std::vector<Subgraph> expected{{0, 1, 2, 3}};

            Finder finder;
            auto actual = find_all_subgraphs(finder, C4);
            expect(name + " recognizes C4", normalize(expected), normalize(actual));
        }

        // 3-0-4
        // | |
        // 2-1-5
        auto G = Graph::from_edges(6, {
            {0, 1}, {1, 2},
            {2, 3}, {3, 0},
            {0, 4}, {1, 5}});

        {
            std::vector<Subgraph> expected({{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 5}, {0, 1, 4, 5}, {0, 2, 3, 4}, {1, 2, 3, 5}});
            Finder finder;
            auto actual = find_all_subgraphs(finder, G);
            expect(name + " recognizes C4P4 in small graph", normalize(expected), normalize(actual));
        }

        {
            Graph forbidden(6);
            forbidden.setEdge({1, 2});
            std::vector<Subgraph> expected({{0, 1, 3, 5}, {0, 1, 4, 5}, {0, 2, 3, 4}});
            Finder finder;
            auto actual = find_all_non_forbidden_subgraphs(finder, G, forbidden);
            expect(name + " recognizes C4P4 in small graph with forbidden edges", normalize(expected), normalize(actual));
        }

        /*    4
              |
            3-0-5
             /|
          b-2-1-6
           /| |\
          a 9 8 7  */
        auto G2 = Graph::from_edges(12, {
            {0, 1}, {1, 2},
            {2, 0}, {0, 3},
            {0, 4}, {0, 5},
            {1, 6}, {1, 7},
            {1, 8}, {2, 9},
            {2, 10}, {2, 11}});

        {
            Finder finder;
            std::vector<Subgraph> expected({{3, 0, 1, 6}, {3, 0, 1, 7}, {3, 0, 1, 8}, {3, 0, 2, 9}, {3, 0, 2, 10}, {3, 0, 2, 11},
                                            {4, 0, 1, 6}, {4, 0, 1, 7}, {4, 0, 1, 8}, {4, 0, 2, 9}, {4, 0, 2, 10}, {4, 0, 2, 11},
                                            {5, 0, 1, 6}, {5, 0, 1, 7}, {5, 0, 1, 8}, {5, 0, 2, 9}, {5, 0, 2, 10}, {5, 0, 2, 11},
                                            {6, 1, 2, 9}, {6, 1, 2, 10}, {6, 1, 2, 11}, {7, 1, 2, 9}, {7, 1, 2, 10}, {7, 1, 2, 11}, {8, 1, 2, 9}, {8, 1, 2, 10}, {8, 1, 2, 11}});
            auto actual = find_all_subgraphs(finder, G2);
            expect(name + " recognizes many P4 in small graph", normalize(expected), normalize(actual));
        }

        {
            auto forbidden = Graph::from_edges(12, {{0, 1}});
            Finder finder;
            std::vector<Subgraph> expected({{3, 0, 2, 9}, {3, 0, 2, 10}, {3, 0, 2, 11},
                                            {4, 0, 2, 9}, {4, 0, 2, 10}, {4, 0, 2, 11},
                                            {5, 0, 2, 9}, {5, 0, 2, 10}, {5, 0, 2, 11},
                                            {6, 1, 2, 9}, {6, 1, 2, 10}, {6, 1, 2, 11}, {7, 1, 2, 9}, {7, 1, 2, 10}, {7, 1, 2, 11}, {8, 1, 2, 9}, {8, 1, 2, 10}, {8, 1, 2, 11}});
            auto actual = find_all_non_forbidden_subgraphs(finder, G2, forbidden);
            expect(name + " recognizes many P4 in small graph with forbidden edges", normalize(expected), normalize(actual));
        }
    }

    template <typename Finder>
    void Finder_finds_P4(const std::string& name) {
        auto G = Graph::make_path_graph(4);

        std::vector<Subgraph> expected{{0, 1, 2, 3}};
        Finder finder;
        auto actual = find_all_subgraphs(finder, G);
        expect(name + " recognizes P4", normalize(expected), normalize(actual));
    }

    template <typename Finder>
    void Finder_finds_P3(const std::string& name) {
        auto G = Graph::make_path_graph(3);

        std::vector<Subgraph> expected{{0, 1, 2}};
        Finder finder;
        auto actual = find_all_subgraphs(finder, G);
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
                auto finder = make_finder(name);

                std::vector<Subgraph> find;
                finder->find(G, [&](Subgraph &&subgraph) {
                    find.push_back(std::move(subgraph));
                    return false;
                });

                std::vector<Subgraph> find_forbidden;
                finder->find(G, F, [&](Subgraph &&subgraph) {
                    find_forbidden.push_back(std::move(subgraph));
                    return false;
                });

                std::vector<std::vector<Subgraph>> find_near;
                for (VertexPair uv : G.vertexPairs()) {
                    std::vector<Subgraph> find_near_uv;
                    finder->find_near(uv, G, [&](Subgraph &&subgraph) {
                        find_near_uv.push_back(std::move(subgraph));
                        return false;
                    });
                    find_near.push_back(normalize(find_near_uv));
                }

                std::vector<std::vector<Subgraph>> find_near_forbidden;
                for (VertexPair uv : G.vertexPairs()) {
                    std::vector<Subgraph> find_near_uv;
                    finder->find_near(uv, G, F, [&](Subgraph &&subgraph) {
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

    template <typename Finder>
    void Finder_finds_C4P4_with_duplicates(const std::string &name) {

        auto to_string = [](auto x) {
            std::stringstream ss;
            ss << x;
            return ss.str();
        };

        {
            auto P4 = Graph::make_path_graph(4);
            Graph forbidden(4);

            std::vector<Subgraph> expected{{0, 1, 2, 3}};
            Finder finder;
            auto actual = find_all_subgraphs_with_duplicates(finder, P4, forbidden);
            expect(name + " recognizes P4 (with duplicates)", expected, actual);
        }

        {
            auto C4 = Graph::make_cycle_graph(4);
            Graph forbidden(4);

            std::vector<Subgraph> expected{{3, 0, 1, 2}, {1, 0, 3, 2}, {0, 1, 2, 3}, {1, 2, 3, 0}};
            Finder finder;
            auto actual = find_all_subgraphs_with_duplicates(finder, C4, forbidden);
            expect(name + " recognizes C4 (with duplicates)", expected, actual);
        }

        {
            auto P4 = Graph::make_path_graph(4);

            VertexPairMap<std::vector<Subgraph>> expected(4);
            expected[{3, 0}] = {{0, 1, 2, 3}};

            for (auto uv : P4.vertexPairs()) {
                auto forbidden = Graph::from_edges(4, {uv});

                Finder finder;
                auto actual = find_all_subgraphs_with_duplicates(finder, P4, forbidden);
                expect(name + " recognizes P4 (with duplicates) forbidden " + to_string(uv), expected[uv], actual);
            }
        }

        {
            auto C4 = Graph::make_cycle_graph(4);

            VertexPairMap<std::vector<Subgraph>> expected(4);
            expected[{0, 1}] = {{1, 2, 3, 0}};  // {{0, 3, 2, 1}}  NOTE(jonas): Consider allowing other orientation
            expected[{1, 2}] = {{1, 0, 3, 2}};  // {{2, 3, 0, 1}}
            expected[{2, 3}] = {{3, 0, 1, 2}};  // {{2, 1, 0, 3}}
            expected[{3, 0}] = {{0, 1, 2, 3}};  // {{3, 2, 1, 0}}

            for (auto uv : C4.vertexPairs()) {
                auto forbidden = Graph::from_edges(4, {uv});

                Finder finder;
                auto actual = find_all_subgraphs_with_duplicates(finder, C4, forbidden);
                expect(name + " recognizes C4 (with duplicates) forbidden " + to_string(uv), expected[uv], actual);
            }
        }
    }

    template <typename Finder>
    void Finder_C4P4_find_near_with_duplicates_1(const std::string &name) {

        auto to_string = [](auto x) {
            std::stringstream ss;
            ss << x;
            return ss.str();
        };

        /*    4
              |
            3-0-5
             /|
          b-2-1-6
           /| |\
          a 9 8 7  */
        auto G2 = Graph::from_edges(12, {
                {0, 1}, {1, 2},
                {2, 0}, {0, 3},
                {0, 4}, {0, 5},
                {1, 6}, {1, 7},
                {1, 8}, {2, 9},
                {2, 10}, {2, 11}});

        VertexPairMap<std::vector<Subgraph>> expected(12);
        std::vector<Subgraph> all_subgraphs = {{3, 0, 1, 6}, {3, 0, 1, 7}, {3, 0, 1, 8}, {3, 0, 2, 9}, {3, 0, 2, 10}, {3, 0, 2, 11},
                                               {4, 0, 1, 6}, {4, 0, 1, 7}, {4, 0, 1, 8}, {4, 0, 2, 9}, {4, 0, 2, 10}, {4, 0, 2, 11},
                                               {5, 0, 1, 6}, {5, 0, 1, 7}, {5, 0, 1, 8}, {5, 0, 2, 9}, {5, 0, 2, 10}, {5, 0, 2, 11},
                                               {6, 1, 2, 9}, {6, 1, 2, 10}, {6, 1, 2, 11}, {7, 1, 2, 9}, {7, 1, 2, 10}, {7, 1, 2, 11}, {8, 1, 2, 9}, {8, 1, 2, 10}, {8, 1, 2, 11}};

        for (auto uv : G2.vertexPairs()) {
            for (Subgraph subgraph : all_subgraphs) {
                if (subgraph.contains(uv)) {
                    auto a = subgraph[0]; auto b = subgraph[1]; auto c = subgraph[2]; auto d = subgraph[3];
                    if (G2.hasEdge({a, d})) {
                        // C4
                        expected[uv].push_back({a, b, c, d});
                        expected[uv].push_back({b, c, d, a});
                        expected[uv].push_back({c, d, a, b});
                        expected[uv].push_back({d, a, b, c});
                    } else {
                        // P4
                        expected[uv].push_back({a, b, c, d});
                    }
                }
            }
        }

        {
            auto forbidden = Graph::from_edges(12, {});
            Finder finder;

            for (auto uv : G2.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near_with_duplicates(uv, G2, forbidden, [&](auto subgraph) {
                    actual.push_back(subgraph);
                    return false;
                });

                expect(name + " find_near_with_duplicates no forbidden " + to_string(uv), normalize(expected[uv]), normalize(actual));
            }
        }
    }

    template <typename Finder>
    void Finder_C4P4_find_near_with_duplicates_2(const std::string &name) {

        auto to_string = [](auto x) {
            std::stringstream ss;
            ss << x;
            return ss.str();
        };

        Finder finder;

        {
            // Output subgraph
            auto G = Graph::make_path_graph(4);
            auto forbidden = Graph::make_empty_graph(4);
            std::vector<Subgraph> expected = {{0, 1, 2, 3}};

            for (auto uv : G.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near_with_duplicates(uv, G, forbidden, [&](Subgraph &&subgraph) {
                    actual.push_back(subgraph);
                    return false;
                });
                expect(name + " find_near_with_duplicates P4 no forbidden " + to_string(uv), expected, actual);
            }
        }

        {
            // Output subgraph, as a-b-c-d has ad = {0, 3} every time.
            auto G = Graph::make_path_graph(4);
            auto forbidden = Graph::from_edges(4, {{0, 3}});

            std::vector<Subgraph> expected = {{0, 1, 2, 3}};

            for (auto uv : G.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near_with_duplicates(uv, G, forbidden, [&](Subgraph &&subgraph) {
                    actual.push_back(subgraph);
                    return false;
                });
                expect(name + " find_near_with_duplicates P4 conversion non-edge forbidden " + to_string(uv), expected, actual);
            }
        }

        {
            // Output no subgraph
            auto G = Graph::make_path_graph(4);
            auto forbidden = Graph::from_edges(4, {{0, 1}});

            std::vector<Subgraph> expected;

            for (auto uv : G.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near_with_duplicates(uv, G, forbidden, [&](Subgraph &&subgraph) {
                    actual.push_back(subgraph);
                    return false;
                });

                expect(name + " find_near_with_duplicates P4 one edge forbidden " + to_string(uv), expected, actual);
            }
        }

        {
            // Output no subgraph
            auto G = Graph::make_path_graph(4);
            auto forbidden = Graph::from_edges(4, {{0, 2}});
            std::vector<Subgraph> expected;

            for (auto uv : G.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near_with_duplicates(uv, G, forbidden, [&](Subgraph &&subgraph) {
                    actual.push_back(subgraph);
                    return false;
                });
                expect(name + " find_near_with_duplicates P4 one non-edge forbidden " + to_string(uv), expected, actual);
            }
        }

        {
            // Output all 4 versions of the subgraph
            auto G = Graph::make_cycle_graph(4);
            auto forbidden = Graph::make_empty_graph(4);
            std::vector<Subgraph> expected = {{0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}};

            for (auto uv : G.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near_with_duplicates(uv, G, forbidden, [&](Subgraph &&subgraph) {
                    subgraph.sortVertices();
                    actual.push_back(subgraph);
                    return false;
                });
                expect(name + " find_near_with_duplicates C4 no forbidden " + to_string(uv), expected, actual);
            }
        }

        {
            // Only output version of a-b-c-d-a with ad = {0, 1}
            auto G = Graph::make_cycle_graph(4);
            auto forbidden = Graph::from_edges(4, {{0, 1}});
            std::vector<Subgraph> expected = {{0, 1, 2, 3}};

            for (auto uv : G.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near_with_duplicates(uv, G, forbidden, [&](Subgraph &&subgraph) {
                    subgraph.sortVertices();
                    actual.push_back(subgraph);
                    return false;
                });
                expect(name + " find_near_with_duplicates C4 one edge forbidden " + to_string(uv), expected, actual);
            }
        }

        {
            auto G = Graph::make_cycle_graph(4);
            auto forbidden = Graph::from_edges(4, {{0, 2}});
            std::vector<Subgraph> expected = {};

            for (auto uv : G.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near_with_duplicates(uv, G, forbidden, [&](Subgraph &&subgraph) {
                    subgraph.sortVertices();
                    actual.push_back(subgraph);
                    return false;
                });
                expect(name + " find_near_with_duplicates C4 one non-edge forbidden " + to_string(uv), expected, actual);
            }
        }
    }

    template <typename Finder>
    void FinderC4P4_for_all_conversionless_edits(const std::string &name) {
        {
            Subgraph subgraph{0, 1, 2, 3};

            std::vector<VertexPair> expected{{0, 1}, {0, 2}, {1, 2}, {1, 3}, {2, 3}};
            Finder finder;

            std::vector<VertexPair> actual;
            finder.for_all_conversionless_edits(subgraph, [&](auto uv) {
                actual.push_back(uv);
                return false;
            });
            expect(name + " lists vertex pairs for C4P4 skipping conversion", expected, actual);
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
        // finders_have_same_output({"NaiveC4P4", "CenterC4P4", "CenterRecC4P4", "EndpointRecC4P4", "NaiveRecC4P4"}, {0, 1});

        Finder_finds_C4P4_with_duplicates<Finder::CenterC4P4>("CenterC4P4");
        FinderC4P4_for_all_conversionless_edits<Finder::CenterC4P4>("CenterC4P4");

        Finder_C4P4_find_near_with_duplicates_1<Finder::CenterC4P4>("CenterC4P4");
        Finder_C4P4_find_near_with_duplicates_2<Finder::CenterC4P4>("CenterC4P4");


        // finders_have_same_output({"CenterRecC5P5", "EndpointRecC5P5", "NaiveRecC5P5"}, {0, 1});

        FinderC4P4_karate_graph_is_solved_with_21_edits();
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FINDERTESTS_H
