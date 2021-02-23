#include <random>

#include "NPS_MWIS_Solver.h"
#include "../../extern/nps_mwis/src/algorithm.h"


namespace nps_mwis {
    // nps_mwis::Solution relies on std::mt19937 generator in the nps_mwis namespace.
    // TODO: Refactor ils_mwis to eliminate global variables
    std::mt19937 generator(0);
}

namespace lower_bound {

    template<Options::FSG SetOfForbiddenSubgraphs>
    Cost NPS_MWIS_Solver<SetOfForbiddenSubgraphs>::calculate_lower_bound(Cost k) {

        // Configuration
        nps_mwis::ArgPack ap;

        ap.complement = false;
        ap.p[0] = 2, ap.p[1] = 4, ap.p[2] = 4, ap.p[3] = 1;
        ap.verbose = false;
        ap.rand_seed = 0;
        ap.target = k + 1;
        ap.iterations = 500;

        auto graph_instance = build_instance();

        if (graph_instance) {

            // nps_mwis does not seem to be able to handle instances with zero vertices.
            if (graph_instance.value().n() == 0)
                return 0;

            // ... or zero edges/one vertex
            if (graph_instance.value().m() == 0) {
                Cost sum = 0;
                for (size_t u = 0; u < graph_instance.value().n(); ++u) {
                    sum += graph_instance.value().weight(u);
                }
                return sum;
            }

            auto solution = nps_mwis::solve(&graph_instance.value(), ap);

            return solution.weight();
        } else {
            return invalid_cost;
        }
    }

    /**
     * Build a Maximum Weight Independent Set (MWIS) instance which encodes the subgraph packing problem.
     *
     * The return type is a ils_mwis Graph.
     *
     * @return
     */
    template<Options::FSG SetOfForbiddenSubgraphs>
    std::optional<nps_mwis::Graph> NPS_MWIS_Solver<SetOfForbiddenSubgraphs>::build_instance() {
        VertexPairMap<std::vector<Vertex>> cliques(m_edit_state->graph().size());

        std::vector<int> weights;

        // The set of subgraphs which share a vertex pair form a clique in the MWIS instance.
        auto exit_state = finder.find(m_edit_state->graph(), [&](Subgraph subgraph) {
            const auto index = weights.size();
            for (auto uv : subgraph.non_converting_edits()) {
                cliques[uv].push_back(index);
            }
            Cost cost = subgraph.calculate_min_cost(m_edit_state->cost_map(), m_edit_state->marked_map());
            weights.push_back(cost);
            // If all vertex pairs are marked then the cost is invalid and the subgraph cannot be destroyed because no
            // vertex pairs can be edited.
            return subgraph_iterators::break_if(cost == invalid_cost);
        });

        if (exit_state == subgraph_iterators::IterationExit::Break)
            return {};

        const auto n = weights.size();

        // The detour through the Graph instance instead directly the ils_mwis::Graph instance is needed because two
        // subgraphs may share multiple vertex pairs and therefore the same edge would be inserted multiple times.
        // ils_mwis::Graph stores a vector of neighbors for each neighbor and does not handle multiple insertions.
        Graph instance_graph(n);

        for (VertexPair uv : m_edit_state->graph().vertex_pairs()) {
            if (m_edit_state->is_marked(uv))
                continue;

            const auto &clique = cliques[uv];
            for (size_t i = 0; i < clique.size(); ++i) {
                for (size_t j = i + 1; j < clique.size(); ++j) {
                    instance_graph.set_edge({clique[i], clique[j]});
                }
            }
        }

        size_t m = 0;
        for (Vertex u : instance_graph.vertices())
            m += instance_graph.degree(u);
        m /= 2;

        // Transform the Graph into a ils_mwis::Graph.
        nps_mwis::Graph ils_mwis_instance_graph(n, m);

        for (size_t u = 0; u < n; ++u) {
            ils_mwis_instance_graph.weight(u) = weights[u];
        }
        for (auto[u, v] : instance_graph.edges()) {
            ils_mwis_instance_graph.addEdge(u, v);
        }

        return ils_mwis_instance_graph;
    }

    template class NPS_MWIS_Solver<Options::FSG::C4P4>;
    template class NPS_MWIS_Solver<Options::FSG::P3>;
}
