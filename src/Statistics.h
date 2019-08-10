//
// Created by jonas on 17.07.19.
//

#ifndef CONCEPT_STATISTICS_H
#define CONCEPT_STATISTICS_H

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#include "definitions.h"

class Statistics {
    using Depth = int;
    Cost min;
    Cost max;
    size_t n_buckets;
    std::vector<int> calls_cost;
    std::vector<int> calls_depth;

public:
    Statistics(Cost min_, Cost max_, size_t n_buckets_) : min(min_), max(max_), n_buckets(n_buckets_), calls_cost(n_buckets + 2), calls_depth(max) {}

    void calls(Depth depth, Cost cost) {
        calls_depth[depth]++;
        calls_cost[bucket(cost)]++;
    }

    [[nodiscard]] std::string yaml() const {
         std::stringstream ss;
         ss << "buckets: [" << min << " " << max << "] " << n_buckets;
         ss << "\ncalls_cost:\n\t";
         for (auto x : calls_cost) ss << x << " ";
         ss << "\ncalls_depth:\n\t";
         for (auto x : calls_depth) ss << x << " ";
         return ss.str();
    }

    friend std::ostream& operator<<(std::ostream &os, const Statistics& statistics) {
        return os << statistics.yaml();
    }

private:
    [[nodiscard]] size_t bucket(Cost cost) const {
        if (cost < min) {
            return 0;
        } else if (cost > max) {
            return n_buckets + 1;
        } else {
            return 1 + n_buckets * (cost - min) / (1 + max - min);
        }
    }
};


#endif //CONCEPT_STATISTICS_H
