#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SOLUTION_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SOLUTION_H


#include "graph/VertexPairMap.h"
#include "graph/Graph.h"
#include "Configuration.h"
#include "forbidden_subgraphs/subgraphs.h"


class Solution {
public:
    Cost cost;
    std::vector<VertexPair> edits;

    Solution(const VertexPairMap<Cost> &costs, std::vector<VertexPair> edits_) : cost(0), edits(std::move(edits_)) {
        for (VertexPair uv : edits) cost += costs[uv];

        std::sort(edits.begin(), edits.end());
    }
    Solution(const Instance &instance, std::vector<VertexPair> edits_) : Solution(instance.costs, std::move(edits_)) {}

    friend std::ostream &operator<<(std::ostream &os, const Solution &solution) {
        os << solution.cost << " { ";
        for (auto uv : solution.edits) os << uv << " ";
        return os << "}";
    }

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const Solution &solution) {
        out << YAML::BeginMap;
        out << YAML::Key << "cost";
        out << YAML::Value << solution.cost;
        out << YAML::Key << "edits";
        out << YAML::Value << YAML::Flow << solution.edits;
        out << YAML::EndMap;
        return out;
    }

    [[nodiscard]] bool is_valid(const Instance &instance, Options::FSG fsg) const {
        Graph graph = instance.graph.copy();

        switch (fsg) {
            case Options::FSG::C4P4: {
                Cost sum{0};
                for (VertexPair uv : edits) {
                    graph.toggle_edge(uv);
                    sum += instance.costs[uv];
                }

                SubgraphT<Options::FSG::C4P4>::Finder finder;
                auto exit_state = finder.find(graph, [&](const auto &) { return subgraph_iterators::IterationControl::Break; });

                return exit_state == subgraph_iterators::IterationExit::Normal;
            }
            default:
                throw std::runtime_error("Solution::is_valid is not specialized for given set of forbidden subgraphs.");
        }
    }

    bool operator<(const Solution &other) const {
        return cost < other.cost;
    }

    bool operator==(const Solution &other) const {
        return cost == other.cost && edits == other.edits;
    }

    bool operator!=(const Solution &other) const {
        return !(*this == other);
    }

    static void filter_inclusion_minimal(std::vector<Solution> &solutions) {
        for (auto &solution : solutions)
            std::sort(solution.edits.begin(), solution.edits.end());

        std::vector<Solution> result;
        for (const auto &a : solutions) {
            bool is_subset = false;
            for (const auto &b : solutions) {
                if (a != b && std::includes(a.edits.begin(), a.edits.end(), b.edits.begin(), b.edits.end())) {
                    is_subset = true;
                    break;
                }
            }
            if (!is_subset)
                result.push_back(a);
        }

        solutions = result;
        std::sort(solutions.begin(), solutions.end(), [](const Solution &a, const Solution &b) {
            if (a.cost == b.cost) {
                return std::lexicographical_compare(a.edits.begin(), a.edits.end(), b.edits.begin(), b.edits.end());
            } else {
                return a.cost < b.cost;
            }
        });
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SOLUTION_H
