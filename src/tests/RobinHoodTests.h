#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ROBINHOODTESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ROBINHOODTESTS_H


#include "robin_hood.h"


class RobinHoodTests {
public:
    RobinHoodTests(int /*seed*/) {}

    void run() {
        std::cout << "\nRobinHoodTests"
                     "\n-----------" << std::endl;

        robin_hood::unordered_map<Subgraph, Cost> map;

        Subgraph a{0, 1, 2, 3};
        map[std::move(a)] = 24;

        for (const auto&[k, v] : map) {
            std::cout << k << " : " << v << "\n";
        }

        std::vector<std::pair<Subgraph, Cost>> subgraphs{{Subgraph{0, 1, 2}, 10}, {Subgraph{3, 4, 5}, 20}, {Subgraph{6, 7, 8}, 30}};
        for (const auto& [subgraph, cost] : subgraphs) {
            std::cout << "(" << subgraph << " , " << cost << ") ";
        }
        std::cout << "\n";
        map.emplace(std::move(subgraphs.back()));
        subgraphs.pop_back();
        map.insert(robin_hood::pair(std::move(subgraphs.back())));
        for (const auto& [subgraph, cost] : subgraphs) {
            std::cout << "(" << subgraph << " , " << cost << ") ";
        }
        std::cout << "\n";

        for (const auto&[k, v] : map) {
            std::cout << k << " : " << v << "\n";
        }
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ROBINHOODTESTS_H
