//
// Created by jonas on 06.02.20.
//

#include <fstream>
#include "LSSWZ_MWIS_SolverLowerBound.h"

Cost LSSWZ_MWIS_SolverLowerBound::calculate_lower_bound(Cost /*k*/) {

    auto instance = build_instance(*finder, m_graph, m_marked, m_costs);
    auto& graph = instance.first;
    auto& vertex_weights = instance.second;

    constexpr auto filename = "lsswz.tmp.graph";

    std::ofstream file(filename);

    size_t n = graph.size();
    size_t m = 0;
    for (auto u : graph.vertices()) {
        m += graph.degree(u);
    }
    m /= 2;

    file << n << " " << m << " " << 10 << "\n";

    for (auto u : graph.vertices()) {
        file << vertex_weights[u] << " ";
        for (auto v : graph.neighbors(u)) {
            file << v + 1 << " ";
        }
        file << "\n";
    }

    file.close();
    exit(0);

    std::string command = "../extern/lsswz_mwis/code/build/weighted_ls lsswz.tmp.graph --time_limit=5 > lsswz.tmp.out";
    auto ret_code = system(command.c_str());
    assert(ret_code >= 0);

    std::ifstream output("lsswz.tmp.out");
    std::string line;
    for (size_t i = 0; i < 5; ++i) {
        std::getline(output, line);
        // std::cout << "line: " << line << "\n";
    }
    std::string _;
    Cost cost = 0;
    std::stringstream ss(line);
    ss >> _ >> cost;


    std::remove(filename);
    std::remove("lsswz.tmp.out");

    std::cout << "cost: " << cost << "\n";

    return cost;
}

std::pair<Graph, std::vector<Cost>>
LSSWZ_MWIS_SolverLowerBound::build_instance(FinderI &finder, const Graph &graph, const VertexPairMap<bool> &marked,
                                            const VertexPairMap<Cost> &costs) {
    VertexPairMap<std::vector<Vertex>> cliques(graph.size());

    std::vector<Cost> weights;

    // The set of subgraphs which share a vertex pair form a clique in the MWIS instance.
    finder.find([&](const Subgraph &subgraph) {
        const auto index = weights.size();
        for (VertexPair uv : subgraph.vertexPairs()) {
            cliques[uv].push_back(index);
        }
        weights.push_back(get_subgraph_cost(subgraph, marked, costs));
        return false;
    });

    const auto n = weights.size();

    Graph instance_graph(n);

    for (VertexPair uv : graph.vertexPairs()) {
        if (marked[uv])
            continue;

        const auto &clique = cliques[uv];
        for (size_t i = 0; i < clique.size(); ++i) {
            for (size_t j = i + 1; j < clique.size(); ++j) {
                instance_graph.setEdge({clique[i], clique[j]});
            }
        }
    }

    return {instance_graph, weights};
}
