//
// Created by jonas on 06.02.20.
//

#include <fstream>
#include "LSSWZ_MWIS_SolverLowerBound.h"

Cost LSSWZ_MWIS_SolverLowerBound::calculate_lower_bound(Cost /*k*/) {

    auto instance = build_instance(*finder, m_graph, m_marked, m_costs);
    if (!instance.has_value()) {
        if (m_config.verbosity > 0)
            std::cout << "lsswz solve    n: - m: - cost: infty  # unsolvable\n";
        return invalid_cost;
    }

    auto& graph = instance->first;
    auto& vertex_weights = instance->second;

    size_t n = graph.size();
    size_t m = 0;
    for (auto u : graph.vertices()) {
        m += graph.degree(u);
    }
    m /= 2;

    // Early exit
    if (graph.size() == 0 || m == 0) {
        Cost cost = 0;
        for (auto u : graph.vertices()) {
            cost += vertex_weights[u];
        }

        if (m_config.verbosity > 0)
            std::cout << "lsswz solve    n: " << graph.size() << " m: " << m << " cost: " << cost << "  # early exit\n";
        return cost;
    }


    // Write instance to file
    constexpr auto instance_filename = "lsswz.tmp.graph";
    std::ofstream file(instance_filename);
    file << n << " " << m << " " << 10 << "\n";

    for (auto u : graph.vertices()) {
        file << vertex_weights[u] << " ";
        for (auto v : graph.neighbors(u)) {
            file << v + 1 << " ";
        }
        file << "\n";
    }
    file.close();


    // Execute MWIS solver
    std::string command = "../extern/lsswz_mwis/code/build/weighted_ls lsswz.tmp.graph"
                          " --console_log --disable_checks --reduction_style dense";
    if (time_limit > 0)
        command += " --time_limit=" + std::to_string(time_limit);
    if (disable_reduction)
        command += " --disable_reduction";

    command += " > lsswz.tmp.out";
    auto ret_code = system(command.c_str());

    if (ret_code != 0) {
        std::remove(instance_filename);
        std::remove("lsswz.tmp.out");
        throw std::runtime_error("weighted_ls exited with a non-zero exit code.");
    }

    // Read cost from output file
    std::ifstream output("lsswz.tmp.out");
    std::string line;
    Cost cost = 0;
    do {
        std::getline(output, line);
        if (m_config.verbosity > 1 && !line.empty()) {
            std::cout << "    " << line << "\n";
        }
        std::string prefix("%MIS_weight ");
        if (!line.compare(0, prefix.size(), prefix)) {
            line = line.substr(prefix.size());
            std::stringstream ss(line);
            ss >> cost;
        }
    } while (!line.empty());

    std::remove(instance_filename);
    std::remove("lsswz.tmp.out");


    /*
     * Am stupid, ignore.
    // The lsswz_mwis seems to sometimes calculate (reduction_offset + 2147483647) instead of (reduction_offset).
    // We try to detect these errors and fix them.
    // TODO: Propose a fix
    constexpr Cost MAX_NodeWeight = std::numeric_limits<int>::max();  // 2147483647; // 2^31 - 1
    static_assert(MAX_NodeWeight > 0, "The MAX_NodeWeight is positive and fits into Cost.");
    if (cost > MAX_NodeWeight / 2) {
        cost -= MAX_NodeWeight;
    }
    */

    if (m_config.verbosity > 0)
        std::cout << "lsswz solve    n: " << graph.size() << " m: " << m << " cost: " << cost << "\n";

    return cost;
}

/**
 * Returns an mwis instance corresponding to the subgraph packing problem.
 *
 * Returns a null option if the subgraph packing problem is not solvable, i.e. all vertex pairs of a forbidden subgraph
 * are marked.
 *
 * @param finder
 * @param graph
 * @param marked
 * @param costs
 * @return
 */
std::optional<std::pair<Graph, std::vector<Cost>>>
LSSWZ_MWIS_SolverLowerBound::build_instance(FinderI &finder, const Graph &graph, const VertexPairMap<bool> &marked,
                                            const VertexPairMap<Cost> &costs) {
    VertexPairMap<std::vector<Vertex>> cliques(graph.size());

    std::vector<Cost> weights;

    // The set of subgraphs which share a vertex pair form a clique in the MWIS instance.
    bool unsolvable = finder.find(graph, [&](const Subgraph &subgraph) {
        const auto index = weights.size();
        for (VertexPair uv : subgraph.vertexPairs()) {
            cliques[uv].push_back(index);
        }
        Cost cost = get_subgraph_cost(subgraph, marked, costs);
        weights.push_back(cost);
        return cost == invalid_cost;
    });

    if (unsolvable)
        return std::nullopt;

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

    return std::make_pair(instance_graph, weights);
}
