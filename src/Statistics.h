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
    using Depth = unsigned int;
    Cost m_min;
    Cost m_max;
    size_t m_num_buckets;

    std::vector<int> m_prunes;
    std::vector<int> m_calls;

public:
    Statistics() : m_min(0), m_max(0), m_num_buckets(0) {}

    Statistics(Cost min, Cost max, size_t num_buckets) : m_min(min), m_max(max), m_num_buckets(num_buckets), m_prunes(m_num_buckets + 2), m_calls(m_num_buckets + 2) {}

    int &calls(Cost cost) {
        return m_calls[bucket(cost)];
    }

    int &prunes(Cost cost) {
        return m_prunes[bucket(cost)];
    }

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const Statistics &statistics) {
        out << YAML::BeginMap;
        out << YAML::Key << "buckets"
            << YAML::Value << YAML::Flow << YAML::BeginSeq;
        out << statistics.m_min << "..." << statistics.m_max;
        out << YAML::EndSeq;
        out << YAML::Key << "calls"
            << YAML::Value << YAML::Flow << statistics.m_calls;
        out << YAML::Key << "prunes"
            << YAML::Value << YAML::Flow << statistics.m_prunes;
        return out << YAML::EndMap;
    }

private:
    [[nodiscard]] size_t bucket(Cost cost) const {
        if (cost < m_min) {
            return 0;
        } else if (cost > m_max) {
            return m_num_buckets + 1;
        } else {
            return 1 + m_num_buckets * (static_cast<unsigned long>(cost - m_min))
                / (static_cast<unsigned long>(1 + m_max - m_min));
        }
    }
};


#endif //CONCEPT_STATISTICS_H
