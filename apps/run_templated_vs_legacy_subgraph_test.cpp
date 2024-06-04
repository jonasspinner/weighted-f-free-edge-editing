#include "../src/forbidden_subgraphs/SubgraphC4P4.h"
#include "../src/graph/GraphIO.h"


#include "../src/tests/FinderTTests.h"


template<Options::FSG Kind>
class Test {
    using Subgraph = SubgraphT<Kind>;
    std::vector<Subgraph> m_vector;
    typename Subgraph::Finder m_finder;
    Graph m_graph;
public:
    Test() : m_graph(Graph::make_cycle_graph(5)) {}

    void fill() {
        m_finder.find(m_graph, [&](Subgraph subgraph) {
            m_vector.push_back(subgraph);
            return subgraph_iterators::IterationControl::Continue;
        });

        for (const auto &x : m_vector) {
            std::cout << x << "\n";
        }

        constexpr auto S = SubgraphT<Options::FSG::C4P4>::C4({0, 1, 2, 3});
        for (auto u : S.vertices()) {
            std::cout << u << " ";
        }
        std::cout << "\n";
        for (auto uv : S.vertex_pairs()) {
            std::cout << uv << " ";
        }
        std::cout << "\n";
        for (auto uv : S.non_converting_edits()) {
            std::cout << uv << " ";
        }
        std::cout << "\n";
    }
};

int main() {
    Test<Options::FSG::C4P4> t;
    t.fill();

    constexpr auto S = SubgraphT<Options::FSG::C4P4>::C4({0, 1, 2, 3});
    std::cout << S << "\n";

    SubgraphT<Options::FSG::C4P4>::Finder finder;
    auto c4 = Graph::make_cycle_graph(4);
    auto c6 = Graph::make_cycle_graph(6);
    auto f1 = Graph::from_edges(6, {{0, 1}});
    auto f2 = Graph::from_edges(6, {{0, 2}});

    std::cout << "---\n";
    finder.find(c4, [](auto subgraph) {
        std::cout << subgraph << "\n";
        return subgraph_iterators::IterationControl::Continue;
    });
    std::cout << "---\n";
    finder.find_unique(c4, [](auto subgraph) {
        std::cout << subgraph << "\n";
        return subgraph_iterators::IterationControl::Continue;
    });
    std::cout << "---\n";
    finder.find(c6, [](auto subgraph) {
        std::cout << subgraph << "\n";
        return subgraph_iterators::IterationControl::Continue;
    });
    std::cout << "---\n";
    finder.find(c6, f1, [](auto subgraph) {
        std::cout << subgraph << "\n";
        return subgraph_iterators::IterationControl::Continue;
    });
    std::cout << "---\n";
    finder.find(c6, f2, [](auto subgraph) {
        std::cout << subgraph << "\n";
        return subgraph_iterators::IterationControl::Continue;
    });
    std::cout << "---\n";
    finder.find_near({0, 1}, c6, f2, [&](auto subgraph) {
        std::cout << subgraph << "\n";
        return subgraph_iterators::IterationControl::Continue;
    });

    std::cout << "---\n";
    auto instance = GraphIO::read_instance("../data/test/bio-nr-4-size-39.graph");
    auto &G = instance.graph;

    std::size_t n{0};
    for (std::size_t i = 0; i < 100; ++i) {
        finder.find(G, [&](auto) {
            ++n;
            return subgraph_iterators::IterationControl::Continue;
        });
    }
    std::cout << n << "\n";

    FinderTTests().run();

    return 0;
}
