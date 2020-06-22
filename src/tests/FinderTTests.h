
#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FINDERTTESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FINDERTTESTS_H


#include "test_utils.h"
#include "../finder/CenterC4P4.h"


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

template<class F>
auto find_unique(F &finder, const Graph& graph) {
    using Subgraph = SubgraphT<Options::FSG::C4P4>;
    std::vector<Subgraph> subgraphs;
    finder.find_unique(graph, [&](Subgraph subgraph) { subgraphs.push_back(subgraph); return false; });
    return subgraphs;
}

template<class F>
auto find(F &finder, const Graph& graph, const Graph& forbidden) {
    using Subgraph = SubgraphT<Options::FSG::C4P4>;
    std::vector<Subgraph> subgraphs;
    finder.find(graph, forbidden, [&](Subgraph subgraph) { subgraphs.push_back(subgraph); return false; });
    return subgraphs;
}

template<class F>
auto find(F &finder, const Graph& graph) {
    using Subgraph = SubgraphT<Options::FSG::C4P4>;
    std::vector<Subgraph> subgraphs;
    finder.find(graph, [&](Subgraph subgraph) { subgraphs.push_back(subgraph); return false; });
    return subgraphs;
}

template<class F>
auto find_unique(F &finder, const Graph& graph, const Graph& forbidden) {
    using Subgraph = SubgraphT<Options::FSG::C4P4>;
    std::vector<Subgraph> subgraphs;
    finder.find_unique(graph, forbidden, [&](Subgraph subgraph) { subgraphs.push_back(subgraph); return false; });
    return subgraphs;
}


class FinderTTests {
    std::mt19937 gen;
public:
    explicit FinderTTests(int seed=0) : gen(static_cast<unsigned long>(seed)) {}


    void Finder_finds_C4() {
        using Subgraph = SubgraphT<Options::FSG::C4P4>;
        using Finder = Subgraph::Finder;

        auto C4 = Graph::make_cycle_graph(4);

        {
            std::vector<Subgraph> expected{Subgraph::C4({1, 0, 3, 2})};

            Finder finder;
            auto actual = find_unique(finder, C4);
            expect("SubgraphT<Options::FSG::C4P4>::Finder recognizes C4", expected, actual);
        }

        // 3-0-4
        // | |
        // 2-1-5
        auto G = Graph::from_edges(6, {
                {0, 1}, {1, 2},
                {2, 3}, {3, 0},
                {0, 4}, {1, 5}});

        {
            std::vector<Subgraph> expected({
                Subgraph::C4({1, 0, 3, 2}),
                Subgraph::P4({4, 0, 1, 2}), Subgraph::P4({3, 0, 1, 5}), Subgraph::P4({4, 0, 1, 5}),
                Subgraph::P4({4, 0, 3, 2}), Subgraph::P4({5, 1, 2, 3})});
            Finder finder;
            auto actual = find_unique(finder, G);
            std::sort(expected.begin(), expected.end());
            std::sort(actual.begin(), actual.end());
            expect("SubgraphT<Options::FSG::C4P4>::Finder recognizes C4P4 in small graph", expected, actual);
        }

        {
            Graph forbidden(6);
            forbidden.setEdge({1, 2});
            std::vector<Subgraph> expected({
                Subgraph::C4({1, 0, 3, 2}),
                Subgraph::P4({3, 0, 1, 5}),
                Subgraph::P4({4, 0, 1, 5}),
                Subgraph::P4({4, 0,  3, 2})});
            Finder finder;
            auto actual = find(finder, G, forbidden);
            std::sort(expected.begin(), expected.end());
            std::sort(actual.begin(), actual.end());
            expect("SubgraphT<Options::FSG::C4P4>::Finder recognizes C4P4 in small graph with forbidden edges", expected, actual);
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
            std::vector<Subgraph> expected({
                Subgraph::P4({3, 0, 1, 6}), Subgraph::P4({3, 0, 1, 7}), Subgraph::P4({3, 0, 1, 8}), Subgraph::P4({3, 0, 2, 9}), Subgraph::P4({3, 0, 2, 10}), Subgraph::P4({3, 0, 2, 11}),
                Subgraph::P4({4, 0, 1, 6}), Subgraph::P4({4, 0, 1, 7}), Subgraph::P4({4, 0, 1, 8}), Subgraph::P4({4, 0, 2, 9}), Subgraph::P4({4, 0, 2, 10}), Subgraph::P4({4, 0, 2, 11}),
                Subgraph::P4({5, 0, 1, 6}), Subgraph::P4({5, 0, 1, 7}), Subgraph::P4({5, 0, 1, 8}), Subgraph::P4({5, 0, 2, 9}), Subgraph::P4({5, 0, 2, 10}), Subgraph::P4({5, 0, 2, 11}),
                Subgraph::P4({6, 1, 2, 9}), Subgraph::P4({6, 1, 2, 10}), Subgraph::P4({6, 1, 2, 11}), Subgraph::P4({7, 1, 2, 9}), Subgraph::P4({7, 1, 2, 10}), Subgraph::P4({7, 1, 2, 11}), Subgraph::P4({8, 1, 2, 9}), Subgraph::P4({8, 1, 2, 10}), Subgraph::P4({8, 1, 2, 11})});
            auto actual = find_unique(finder, G2);
            std::sort(expected.begin(), expected.end());
            std::sort(actual.begin(), actual.end());
            expect("SubgraphT<Options::FSG::C4P4>::Finder recognizes many P4 in small graph", expected, actual);
        }

        {
            auto forbidden = Graph::from_edges(12, {{0, 1}});
            Finder finder;
            std::vector<Subgraph> expected({
                Subgraph::P4({3, 0, 2, 9}), Subgraph::P4({3, 0, 2, 10}), Subgraph::P4({3, 0, 2, 11}),
                Subgraph::P4({4, 0, 2, 9}), Subgraph::P4({4, 0, 2, 10}), Subgraph::P4({4, 0, 2, 11}),
                Subgraph::P4({5, 0, 2, 9}), Subgraph::P4({5, 0, 2, 10}), Subgraph::P4({5, 0, 2, 11}),
                Subgraph::P4({6, 1, 2, 9}), Subgraph::P4({6, 1, 2, 10}), Subgraph::P4({6, 1, 2, 11}), Subgraph::P4({7, 1, 2, 9}), Subgraph::P4({7, 1, 2, 10}), Subgraph::P4({7, 1, 2, 11}), Subgraph::P4({8, 1, 2, 9}), Subgraph::P4({8, 1, 2, 10}), Subgraph::P4({8, 1, 2, 11})});
            auto actual = find(finder, G2, forbidden);
            std::sort(expected.begin(), expected.end());
            std::sort(actual.begin(), actual.end());
            expect("SubgraphT<Options::FSG::C4P4>::Finder recognizes many P4 in small graph with forbidden edges", expected, actual);
        }
    }


    void Finder_finds_P4() {
        using Subgraph = SubgraphT<Options::FSG::C4P4>;
        using Finder = Subgraph::Finder;

        auto G = Graph::make_path_graph(4);

        std::vector<Subgraph> expected{Subgraph::P4({0, 1, 2, 3})};
        Finder finder;
        auto actual = find(finder, G);
        expect("SubgraphT<Options::FSG::C4P4>::Finder recognizes P4", expected, actual);
    }


    void Finder_is_seed_independent(const std::vector<int> &seeds) {
        using Subgraph = SubgraphT<Options::FSG::C4P4>;
        using Finder = Subgraph::Finder;

        Graph G = random_graph(10, 40, gen);

        auto n = [](Subgraph subgraph) {
            auto r = subgraph.vertices();
            std::vector<Vertex> v(r.begin(), r.end());
            std::sort(v.begin(), v.end());
            return Subgraph::P4({v[0], v[1], v[2], v[3]});
        };

        std::vector<std::pair<int, std::vector<Subgraph>>> results;
        for (auto seed : seeds) {
            auto P = Permutation(G.size(), seed);
            auto P_r = P.reverse();

            auto G_p = P[G];

            std::vector<Subgraph> subgraphs;
            Finder finder;
            finder.find(G_p, [&](Subgraph subgraph) {
                subgraphs.push_back(n(subgraph));
                return false;
            });

            results.emplace_back(seed, std::move(subgraphs));
        }

        for (size_t i = 0; i < seeds.size(); ++i) {
            for (size_t j = i + 1; j < seeds.size(); ++j) {
                const auto &[seed_i, subgraphs_i] = results[i];
                const auto &[seed_j, subgraphs_j] = results[j];
                std::stringstream ss;
                ss << "Finder SubgraphT<Options::FSG::C4P4>::Finder finds the same subgraphs with seed " << seed_i << " and " << seed_j;
                expect(ss.str(), subgraphs_i, subgraphs_j);
            }
        }
    }


    void Finder_finds_C4P4_with_duplicates() {
        using Subgraph = SubgraphT<Options::FSG::C4P4>;
        using Finder = Subgraph::Finder;

        auto to_string = [](auto x) {
            std::stringstream ss;
            ss << x;
            return ss.str();
        };

        {
            auto P4 = Graph::make_path_graph(4);
            Graph forbidden(4);

            std::vector<Subgraph> expected{Subgraph::P4({0, 1, 2, 3})};
            Finder finder;
            auto actual = find(finder, P4, forbidden);
            expect("SubgraphT<Options::FSG::C4P4>::Finder recognizes P4 (with duplicates)", expected, actual);
        }

        {
            auto C4 = Graph::make_cycle_graph(4);
            Graph forbidden(4);

            std::vector<Subgraph> expected{Subgraph::C4({3, 0, 1, 2}), Subgraph::C4({1, 0, 3, 2}), Subgraph::C4({0, 1, 2, 3}), Subgraph::C4({1, 2, 3, 0})};
            Finder finder;
            auto actual = find(finder, C4, forbidden);
            expect("SubgraphT<Options::FSG::C4P4>::Finder recognizes C4 (with duplicates)", expected, actual);
        }

        {
            auto P4 = Graph::make_path_graph(4);

            VertexPairMap<std::vector<Subgraph>> expected(4);
            expected[{3, 0}] = {Subgraph::P4({0, 1, 2, 3})};

            for (auto uv : P4.vertexPairs()) {
                auto forbidden = Graph::from_edges(4, {uv});

                Finder finder;
                auto actual = find(finder, P4, forbidden);
                expect("SubgraphT<Options::FSG::C4P4>::Finder recognizes P4 (with duplicates) forbidden " + to_string(uv), expected[uv], actual);
            }
        }

        {
            auto C4 = Graph::make_cycle_graph(4);

            VertexPairMap<std::vector<Subgraph>> expected(4);
            expected[{0, 1}] = {Subgraph::C4({1, 2, 3, 0})};  // {{0, 3, 2, 1}}  NOTE(jonas): Consider allowing other orientation
            expected[{1, 2}] = {Subgraph::C4({1, 0, 3, 2})};  // {{2, 3, 0, 1}}
            expected[{2, 3}] = {Subgraph::C4({3, 0, 1, 2})};  // {{2, 1, 0, 3}}
            expected[{3, 0}] = {Subgraph::C4({0, 1, 2, 3})};  // {{3, 2, 1, 0}}

            for (auto uv : C4.vertexPairs()) {
                auto forbidden = Graph::from_edges(4, {uv});

                Finder finder;
                auto actual = find(finder, C4, forbidden);
                expect("SubgraphT<Options::FSG::C4P4>::Finder recognizes C4 (with duplicates) forbidden " + to_string(uv), expected[uv], actual);
            }
        }
    }


    void Finder_C4P4_find_near_with_duplicates_1() {
        using Subgraph = SubgraphT<Options::FSG::C4P4>;
        using Finder = Subgraph::Finder;

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
        std::vector<Subgraph> all_subgraphs = {
           Subgraph::P4({3, 0, 1, 6}), Subgraph::P4({3, 0, 1, 7}), Subgraph::P4({3, 0, 1, 8}), Subgraph::P4({3, 0, 2, 9}), Subgraph::P4({3, 0, 2, 10}), Subgraph::P4({3, 0, 2, 11}),
           Subgraph::P4({4, 0, 1, 6}), Subgraph::P4({4, 0, 1, 7}), Subgraph::P4({4, 0, 1, 8}), Subgraph::P4({4, 0, 2, 9}), Subgraph::P4({4, 0, 2, 10}), Subgraph::P4({4, 0, 2, 11}),
           Subgraph::P4({5, 0, 1, 6}), Subgraph::P4({5, 0, 1, 7}), Subgraph::P4({5, 0, 1, 8}), Subgraph::P4({5, 0, 2, 9}), Subgraph::P4({5, 0, 2, 10}), Subgraph::P4({5, 0, 2, 11}),
           Subgraph::P4({6, 1, 2, 9}), Subgraph::P4({6, 1, 2, 10}), Subgraph::P4({6, 1, 2, 11}), Subgraph::P4({7, 1, 2, 9}), Subgraph::P4({7, 1, 2, 10}), Subgraph::P4({7, 1, 2, 11}), Subgraph::P4({8, 1, 2, 9}), Subgraph::P4({8, 1, 2, 10}), Subgraph::P4({8, 1, 2, 11})};

        for (auto uv : G2.vertexPairs()) {
            for (Subgraph subgraph : all_subgraphs) {
                if (subgraph.contains(uv)) {
                    auto a = subgraph[0]; auto b = subgraph[1]; auto c = subgraph[2]; auto d = subgraph[3];
                    if (G2.hasEdge({a, d})) {
                        // C4
                        expected[uv].push_back(Subgraph::C4({a, b, c, d}));
                        expected[uv].push_back(Subgraph::C4({b, c, d, a}));
                        expected[uv].push_back(Subgraph::C4({c, d, a, b}));
                        expected[uv].push_back(Subgraph::C4({d, a, b, c}));
                    } else {
                        // P4
                        expected[uv].push_back(Subgraph::P4({a, b, c, d}));
                    }
                }
            }
        }

        {
            auto forbidden = Graph::from_edges(12, {});
            Finder finder;

            for (auto uv : G2.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near(uv, G2, forbidden, [&](auto subgraph) {
                    actual.push_back(subgraph);
                    return false;
                });

                expect("SubgraphT<Options::FSG::C4P4>::Finder find_near no forbidden " + to_string(uv), expected[uv], actual);
            }
        }
    }

    void Finder_C4P4_find_near_with_duplicates_2() {
        using Subgraph = SubgraphT<Options::FSG::C4P4>;
        using Finder = Subgraph::Finder;

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
            std::vector<Subgraph> expected = {Subgraph::P4({0, 1, 2, 3})};

            for (auto uv : G.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near(uv, G, forbidden, [&](Subgraph subgraph) {
                    actual.push_back(subgraph);
                    return false;
                });
                expect("SubgraphT<Options::FSG::C4P4>::Finder find_near P4 no forbidden " + to_string(uv), expected, actual);
            }
        }

        {
            // Output subgraph, as a-b-c-d has ad = {0, 3} every time.
            auto G = Graph::make_path_graph(4);
            auto forbidden = Graph::from_edges(4, {{0, 3}});

            std::vector<Subgraph> expected = {Subgraph::P4({0, 1, 2, 3})};

            for (auto uv : G.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near(uv, G, forbidden, [&](Subgraph subgraph) {
                    actual.push_back(subgraph);
                    return false;
                });
                expect("SubgraphT<Options::FSG::C4P4>::Finder find_near P4 conversion non-edge forbidden " + to_string(uv), expected, actual);
            }
        }

        {
            // Output no subgraph
            auto G = Graph::make_path_graph(4);
            auto forbidden = Graph::from_edges(4, {{0, 1}});

            std::vector<Subgraph> expected;

            for (auto uv : G.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near(uv, G, forbidden, [&](Subgraph subgraph) {
                    actual.push_back(subgraph);
                    return false;
                });

                expect("SubgraphT<Options::FSG::C4P4>::Finder find_near P4 one edge forbidden " + to_string(uv), expected, actual);
            }
        }

        {
            // Output no subgraph
            auto G = Graph::make_path_graph(4);
            auto forbidden = Graph::from_edges(4, {{0, 2}});
            std::vector<Subgraph> expected;

            for (auto uv : G.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near(uv, G, forbidden, [&](Subgraph subgraph) {
                    actual.push_back(subgraph);
                    return false;
                });
                expect("SubgraphT<Options::FSG::C4P4>::Finder find_near P4 one non-edge forbidden " + to_string(uv), expected, actual);
            }
        }

        {
            // Output all 4 versions of the subgraph
            auto G = Graph::make_cycle_graph(4);
            auto forbidden = Graph::make_empty_graph(4);
            std::vector<Subgraph> expected = {
                Subgraph::C4({0, 1, 2, 3}),
                Subgraph::C4({1, 2, 3, 0}),
                Subgraph::C4({ 1, 0, 3, 2}),
                Subgraph::C4({3, 0, 1, 2})};

            for (auto uv : G.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near(uv, G, forbidden, [&](Subgraph subgraph) {
                    actual.push_back(subgraph);
                    return false;
                });
                std::sort(expected.begin(), expected.end());
                std::sort(actual.begin(), actual.end());
                expect("SubgraphT<Options::FSG::C4P4>::Finder find_near C4 no forbidden " + to_string(uv), expected, actual);
            }
        }

        {
            // Only output version of a-b-c-d-a with ad = {0, 1}
            auto G = Graph::make_cycle_graph(4);
            auto forbidden = Graph::from_edges(4, {{0, 1}});
            std::vector<Subgraph> expected = {Subgraph::C4({1, 2, 3, 0})};

            for (auto uv : G.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near(uv, G, forbidden, [&](Subgraph subgraph) {
                    actual.push_back(subgraph);
                    return false;
                });
                expect("SubgraphT<Options::FSG::C4P4>::Finder find_near C4 one edge forbidden " + to_string(uv), expected, actual);
            }
        }

        {
            auto G = Graph::make_cycle_graph(4);
            auto forbidden = Graph::from_edges(4, {{0, 2}});
            std::vector<Subgraph> expected = {};

            for (auto uv : G.vertexPairs()) {
                std::vector<Subgraph> actual;
                finder.find_near(uv, G, forbidden, [&](Subgraph subgraph) {
                    actual.push_back(subgraph);
                    return false;
                });
                expect("SubgraphT<Options::FSG::C4P4>::Finder find_near C4 one non-edge forbidden " + to_string(uv), expected, actual);
            }
        }
    }

    void FinderC4P4_for_all_conversionless_edits() {
        using Subgraph = SubgraphT<Options::FSG::C4P4>;
        using Finder = Subgraph::Finder;

        {
            Subgraph subgraph = Subgraph::C4({0, 1, 2, 3});

            std::vector<VertexPair> expected{{0, 1}, {0, 2}, {1, 2}, {1, 3}, {2, 3}};
            Finder finder;

            std::vector<VertexPair> actual;
            for (auto uv : subgraph.non_converting_edits()) {
                actual.push_back(uv);
            }
            expect("SubgraphT<Options::FSG::C4P4>::Finder lists vertex pairs for C4P4 skipping conversion", expected, actual);
        }
    }
    
    void Finder_has_same_ouput_as_VectorSubgraph() {
        using Subgraph = SubgraphT<Options::FSG::C4P4>;

        auto ensure_direction = [](auto &vertices) {
            if (vertices[1] > vertices[2]) {
                std::swap(vertices[0], vertices[3]);
                std::swap(vertices[1], vertices[2]);
            }
        };

        auto ensure_rotation = [](auto &vertices) {
            auto [a, b, c, d] = vertices;
            if (a < std::min({b, c, d})) {
                vertices = {d, a, b, c};
            } else if (b < std::min({a, c, d})) {
                vertices = {a, b, c, d};
            } else if (c < std::min({a, b, d})) {
                vertices = {b, c, d, a};
            } else {
                vertices = {c, d, a, b};
            }
            if (vertices[0] > vertices[2]) {
                std::swap(vertices[0], vertices[2]);
            }
        };

        auto instance = GraphIO::read_instance("../data/bio/bio-nr-4-size-39.graph");
        auto &G = instance.graph;

        auto E = Graph::make_empty_graph(G.size());

        auto F = Graph::make_empty_graph(G.size());
        std::uniform_int_distribution<Vertex> dist(0, G.size() - 1);
        for (std::size_t i = 0; i < 2 * G.size(); ++i) {
            auto u = dist(gen), v = dist(gen);
            if (u == v) continue;
            F.setEdge({u, v});
        }

        Subgraph::Finder finder1;
        ::Finder::CenterC4P4 finder2;

        {
            std::vector<Subgraph> output1, output2;

            finder1.find(G, [&](auto subgraph) {
                output1.push_back(subgraph);
                return false;
            });

            finder2.find_with_duplicates(G, [&](auto vector_subgraph) {
                std::array<Vertex, 4> v{vector_subgraph[0], vector_subgraph[1], vector_subgraph[2], vector_subgraph[3]};
                ensure_direction(v);
                if (G.hasEdge({v[0], v[3]})) {
                    output2.push_back(Subgraph::C4(v));
                } else {
                    output2.push_back(Subgraph::P4(v));
                }
                return false;
            });

            std::sort(output1.begin(), output1.end());
            std::sort(output2.begin(), output2.end());
            expect("SubgraphT<Options::FSG::C4P4>::Finder and ::Finder::CenterC4P4 have same find output", output1, output2);
        }

        {
            std::vector<Subgraph> output1, output2;

            finder1.find(G, F,[&](auto subgraph) {
                output1.push_back(subgraph);
                return false;
            });

            finder2.find_with_duplicates(G, F, [&](auto vector_subgraph) {
                std::array<Vertex, 4> v{vector_subgraph[0], vector_subgraph[1], vector_subgraph[2], vector_subgraph[3]};
                ensure_direction(v);
                if (G.hasEdge({v[0], v[3]})) {
                    output2.push_back(Subgraph::C4(v));
                } else {
                    output2.push_back(Subgraph::P4(v));
                }
                return false;
            });

            std::sort(output1.begin(), output1.end());
            std::sort(output2.begin(), output2.end());
            expect("SubgraphT<Options::FSG::C4P4>::Finder and ::Finder::CenterC4P4 have same find output with forbidden pairs", output1, output2);
        }

        {
            std::vector<Subgraph> output1, output2;

            finder1.find_unique(G, [&](auto subgraph) {
                output1.push_back(subgraph);
                return false;
            });

            finder2.find(G, [&](auto vector_subgraph) {
                std::array<Vertex, 4> v{vector_subgraph[0], vector_subgraph[1], vector_subgraph[2], vector_subgraph[3]};
                ensure_direction(v);
                if (G.hasEdge({v[0], v[3]})) {
                    ensure_rotation(v);
                    assert(Subgraph::is_valid_C4(G, E, v));
                    output2.push_back(Subgraph::C4(v));
                } else {
                    assert(Subgraph::is_valid_P4(G, E, v));
                    output2.push_back(Subgraph::P4(v));
                }
                return false;
            });

            std::sort(output1.begin(), output1.end());
            std::sort(output2.begin(), output2.end());
            expect("SubgraphT<Options::FSG::C4P4>::Finder and ::Finder::CenterC4P4 have same find_unique output", output1, output2);
        }

        {
            std::vector<std::vector<Subgraph>> output1, output2;

            for (auto uv : G.vertexPairs()) {
                output1.emplace_back();
                finder1.find_near(uv, G, F, [&](auto subgraph) {
                    output1.back().push_back(subgraph);
                    return false;
                });

                output2.emplace_back();
                finder2.find_near_with_duplicates(uv, G, F, [&](auto vector_subgraph) {
                    std::array<Vertex, 4> v{vector_subgraph[0], vector_subgraph[1], vector_subgraph[2], vector_subgraph[3]};
                    ensure_direction(v);
                    if (G.hasEdge({v[0], v[3]})) {
                        output2.back().push_back(Subgraph::C4(v));
                    } else {
                        output2.back().push_back(Subgraph::P4(v));
                    }
                    return false;
                });

                std::sort(output1.back().begin(), output1.back().end());
                std::sort(output2.back().begin(), output2.back().end());
            }
            expect("SubgraphT<Options::FSG::C4P4>::Finder and ::Finder::CenterC4P4 have same find_near output with forbidden vertex pairs", output1, output2);
        }
    }

    void run() {
        std::cout << "\nFinderTTests"
                     "\n-----------" << std::endl;


        Finder_finds_C4();

        Finder_finds_C4P4_with_duplicates();
        FinderC4P4_for_all_conversionless_edits();

        Finder_C4P4_find_near_with_duplicates_1();
        Finder_C4P4_find_near_with_duplicates_2();

        Finder_has_same_ouput_as_VectorSubgraph();
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FINDERTTESTS_H


