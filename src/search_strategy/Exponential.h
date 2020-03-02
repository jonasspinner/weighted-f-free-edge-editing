//
// Created by jonas on 27.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_EXPONENTIAL_H
#define WEIGHTED_F_FREE_EDGE_EDITING_EXPONENTIAL_H


#include <deque>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "SearchStrategyI.h"


/**
 * From initial implementation as method:
 *
 * This search strategy makes the assumption that the number of calls grows exponentially with the parameter k.
 * We want to increment the paramter k in such a way that each step the number of calls doubles. If we can achieve
 * that, then the overall amount of work (=number of calls) is dominated by the last step.
 *
 *  We try to estimate the relationship between k and the number of calls and fit a log linear model with the
 *  information of the last calls (max_history_length).
 *
 *  In early steps with not enough information the result of the model is very unstable. We fix that by clamping the
 *  size of the step between search_delta estimate (quantile) and a multiple of the last step size (max_step_factor).
 *
 *  TODO: Add clamping.
 */
class Exponential : public SearchStrategyI {
    std::deque<Cost> m_cost_history;
    std::deque<size_t> m_num_calls_history;

    size_t m_max_history_size;
    double m_targeted_growth_factor;

public:
    void set_initial_search_k(Cost initial_search_k) override {
        push_history(initial_search_k);
    }

    std::optional<Cost> get_next_search_k() override {
        Cost next_search_k = find_next_k(m_cost_history, m_num_calls_history, m_targeted_growth_factor);

        push_history(next_search_k);

        return next_search_k;
    }

    void register_call(Cost /*k*/) override {
        assert(!m_num_calls_history.empty());
        ++m_num_calls_history.back();
    }

    void register_bound(Cost /*k*/, Cost /*lower_bound_k*/) override {}

private:
    void push_history(Cost k) {
        if (m_cost_history.size() == m_max_history_size) {
            m_cost_history.pop_front();
            m_num_calls_history.pop_front();
        }

        m_cost_history.push_back(k);
        m_num_calls_history.push_back(0);
    }

    static Cost find_next_k(const std::deque<Cost> &ks, const std::deque<size_t> &calls, double factor) {

        // Convert ks to double
        std::vector<double> xs(ks.begin(), ks.end());

        // Convert calls to double and apply log-transform
        std::vector<double> ys;
        std::transform(calls.begin(), calls.end(), std::back_inserter(ys),
                       [](size_t c) { return std::log(static_cast<double>(c)); });

        const auto[alpha, beta] = solve_linear_regression(xs, ys);

        double next_num_calls = std::log(factor * calls.back());
        double next_x = (next_num_calls - alpha) / beta;

        return std::ceil(next_x);
    }

    /**
     * Solves linear regression model of the form
     *
     *      y ~ alpha + beta * x
     *
     *  for alpha and beta.
     *
     *  References:
     *      [1] https://en.wikipedia.org/wiki/Simple_linear_regression
     *
     * @param xs
     * @param ys
     * @return
     */
    static std::pair<double, double>
    solve_linear_regression(const std::vector<double> &xs, const std::vector<double> &ys) {
        double x_sum = 0, y_sum = 0;
        for (size_t i = 0; i < xs.size(); ++i) {
            x_sum += xs[i];
            y_sum += ys[i];
        }
        auto x_mean = x_sum / xs.size();
        auto y_mean = y_sum / ys.size();

        double s_xy = 0, s2_x = 0;
        for (size_t i = 0; i < xs.size(); ++i) {
            s_xy += (xs[i] - x_mean) * (ys[i] - y_mean);
            s2_x += (xs[i] - x_mean) * (xs[i] - x_mean);
        }

        double beta = s_xy / s2_x;
        double alpha = y_mean - beta * x_mean;

        return {alpha, beta};
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_EXPONENTIAL_H
