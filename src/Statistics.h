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
    std::vector<Cost> m_seperators;

    std::vector<int> m_prunes;
    std::vector<int> m_calls;

public:
    Statistics() : m_min(0), m_max(0), m_num_buckets(0) {}

    /**
     * A class for keeping statistics on the amount of calls and prunes per remaining editing cost.
     *
     * The values are put into 'num_buckets' buckets spanning (-\infty, max]. 'num_buckets - 1' are evenly spaced in
     * (min, max] and an extra bucket for (-\infty, min] is added.
     *
     * @param min
     * @param max
     * @param num_buckets
     */
    Statistics(Cost min, Cost max, size_t num_buckets) : m_min(min), m_max(max),
            m_num_buckets(std::min(num_buckets, static_cast<size_t>(m_max - m_min))), m_seperators(m_num_buckets),
            m_prunes(m_num_buckets), m_calls(m_num_buckets) {

        for (size_t i = 0; i < m_seperators.size(); ++i) {
            m_seperators[i] = m_min + static_cast<Cost>(i) * (m_max - m_min) / static_cast<Cost>(m_num_buckets - 1);
        }

        if (m_num_buckets == 0) {
            m_seperators.push_back(m_max);
            m_prunes.push_back(0);
            m_calls.push_back(0);
        }

        /*
        for (size_t i = 0; i < m_seperators.size() - 1; ++i) {
            auto x = m_seperators[i];
            std::cout << x - 1 << " " << bucket(x - 1) << "\n";
            std::cout << x << " " << bucket(x) << "\n";
            std::cout << x + 1 << " " << bucket(x + 1) << "\n";
        }
        std::cout << m_seperators.back() - 1 << " " << bucket(m_seperators.back() - 1) << "\n";
        std::cout << m_seperators.back() << " " << bucket(m_seperators.back()) << "\n";
        */
    }

    int &calls(Cost cost) {
        assert(!m_calls.empty());
        return m_calls[bucket(cost)];
    }

    int &prunes(Cost cost) {
        assert(!m_prunes.empty());
        return m_prunes[bucket(cost)];
    }

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const Statistics &statistics) {
        using namespace YAML;
        out << BeginMap;
        out << Key << "min" << Value << statistics.m_min;
        out << Key << "max" << Value << statistics.m_max;
        out << Key << "num_buckets" << Value << statistics.m_num_buckets;
        out << Key << "seperators" << Value << Flow << statistics.m_seperators;

        out << Key << "buckets" << Value << Flow << BeginSeq;
        std::stringstream ss;
        const auto& seps = statistics.m_seperators;
        ss << "(-\\infty, " << seps[0] << "]";
        out << ss.str();
        for (size_t i = 0; i < seps.size() - 1; ++i) {
            ss.clear();
            ss.str(std::string());
            ss << "(" << seps[i] << ", " << seps[i + 1] << "]";
            out << ss.str();
        }
        out << EndSeq;

        out << Key << "calls"
            << Value << Flow << statistics.m_calls;
        out << Key << "prunes"
            << Value << Flow << statistics.m_prunes;
        return out << EndMap;
    }

    friend std::ostream &operator<<(std::ostream &os, const Statistics &statistics) {
        YAML::Emitter out;
        out << statistics;
        return os << out.c_str();
    }

private:
    [[nodiscard]] size_t bucket(Cost cost) const {
        assert(cost <= m_max);
        if (cost <= m_min) {
            return 0;
        } else if (cost > m_max) {
            throw std::runtime_error("invalid cost");
        } else {
            return 1 + (m_num_buckets - 1) * (static_cast<unsigned long>(cost - m_min))
                / (static_cast<unsigned long>(1 + m_max - m_min));
        }
    }
};


#endif //CONCEPT_STATISTICS_H
