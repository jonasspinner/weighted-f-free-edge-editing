//
// Created by jonas on 19.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SOLUTION_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SOLUTION_H


#include "graph/VertexPairMap.h"
#include "graph/Graph.h"
#include "Configuration.h"
#include "interfaces/FinderI.h"
#include "finder/util.h"


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
        out << YAML::Value << solution.edits;
        out << YAML::EndMap;
        return out;
    }

    [[nodiscard]] bool is_valid(const Instance &instance, Options::FSG fsg) const {
        Graph graph(instance.graph);
        std::unique_ptr<FinderI> finder = Finder::make(fsg, graph);
        Cost sum = 0;

        for (VertexPair uv : edits) {
            graph.toggle_edge(uv);
            sum += instance.costs[uv];
        }

        bool found_forbidden_subgraph = finder->find([&](const Subgraph &) { return true; });

        return !found_forbidden_subgraph;
    }

    bool operator<(const Solution &other) const {
        return cost < other.cost;
    }

    bool operator==(const Solution &other) const {
        return cost == other.cost && edits == other.edits;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_SOLUTION_H
