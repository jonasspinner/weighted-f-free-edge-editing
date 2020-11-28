#ifndef WEIGHTED_F_FREE_EDGE_EDITING_FPTSOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_FPTSOLVER_H


#include <deque>

#include "Solver.h"
#include "../Configuration.h"
#include "../editor/Editor.h"

class FPTSolver : public Solver {
public:
    struct Stat {
        Cost k;
        int calls;
        double time;
        Stat(Cost k_, int calls_, double time_) : k(k_), calls(calls_), time(time_) {}
    };

private:
    Configuration m_config;
    std::vector<Stat> m_stats;
public:
    explicit FPTSolver(Configuration config) : m_config(std::move(config)) {}

    /**
     * Solves the given instance.
     *
     * @param instance
     * @return
     */
    Result solve(const Instance &instance) override {
        switch (m_config.search_strategy) {
            case Options::FPTSearchStrategy::PrunedDelta:
            case Options::FPTSearchStrategy::Exponential:
            case Options::FPTSearchStrategy::IncrementByMinCost:
            case Options::FPTSearchStrategy::IncrementByMultiplier:
                if (!m_config.find_all_solutions) {
                    std::cerr << "--search-strategy " << m_config.search_strategy << " needs --all 1" << std::endl;
                    abort();
                }
                break;
            case Options::FPTSearchStrategy::Fixed:
            default:
                break;
        }
        switch (m_config.search_strategy) {
            case Options::FPTSearchStrategy::Fixed:
                return edit_fixed(m_config.k_max, instance, m_config, m_stats);
            case Options::FPTSearchStrategy::PrunedDelta:
                return search_delta(instance, m_config, m_stats);
            case Options::FPTSearchStrategy::Exponential:
                return search_exponential(instance, m_config, m_stats);
            case Options::FPTSearchStrategy::IncrementByMinCost:
            {
                Cost delta = std::numeric_limits<Cost>::max();
                for (VertexPair uv : instance.graph.vertexPairs()) {
                    if (0 < instance.costs[uv] && instance.costs[uv] < delta) {
                        delta = instance.costs[uv];
                    }
                }
                return search_incremental(instance, m_config, m_stats, delta);
            }
            case Options::FPTSearchStrategy::IncrementByMultiplier:
                return search_incremental(instance, m_config, m_stats, std::ceil(m_config.multiplier));
            default:
                return Result::Unsolved();
        }
    }

    [[nodiscard]] const std::vector<Stat> &stats() const {
        return m_stats;
    }

private:
    /**
     * Search for solutions L with c(L) <= k.
     *
     * @param k
     * @param instance
     * @param config
     * @param stats
     * @return
     */
    static Result edit_fixed(Cost k, const Instance &instance, const Configuration &config, std::vector<Stat> &stats) {
        if (k < 0)
            return Result::Unsolved();
        std::vector<Solution> solutions;

        Editor editor(instance.graph.copy(), instance.costs, config);

        auto k_min = editor.initial_lower_bound();
        if (k_min > k)
            return Result::Unsolved();

        auto start = std::chrono::steady_clock::now();

        bool solved = editor.edit(k, [&](const std::vector<VertexPair> &edits) {
            solutions.emplace_back(instance, edits);
        }, [](Cost, Cost) {});

        auto end = std::chrono::steady_clock::now();

        stats.emplace_back(k, editor.stats().allCalls(), std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        if (solved) {
            std::sort(solutions.begin(), solutions.end());
            Solution::filter_inclusion_minimal(solutions);
            return Result::Solved(solutions);
        } else {
            return Result::Unsolved();
        }
    }

    /**
     * Searches for a solution with minimal solution cost.
     *
     * For each edit step we record the following:
     * + when a branch is pruned, what could increment in available editing
     *   cost would be necessary that the branch would not be pruned for the calculated lower bound?
     *   We call this value delta.
     * + From all delta choose the value from the specified quantile.
     * + Increment k by the chosen value.
     *
     * Caveat: Lower bounds on the editing cost are no longer maximized, when they have already reached the given value.
     *
     *
     * @param instance
     * @param config
     * @param quantile The quantile of delta values.
     * @return
     */
    static Result search_delta(const Instance &instance, const Configuration &config,
            std::vector<Stat> &stats, double quantile = 0.5) {
        Editor editor(instance.graph.copy(), instance.costs, config);

        std::chrono::seconds timelimit(config.timelimit);

        Cost k = editor.initial_lower_bound();
        bool solved;
        std::vector<Solution> solutions;

        do {
            Cost min_delta = std::numeric_limits<Cost>::max();
            std::vector<Cost> deltas;
            size_t num_calls = 0;
            Cost k_min = 0;

            auto result_cb = [&](const std::vector<VertexPair> &edits) {
                solutions.emplace_back(instance, edits);
                k_min = k - solutions.back().cost;
            };

            auto prune_cb = [&](Cost pruned_k, Cost lb) {
                if (lb == invalid_cost) return;
                min_delta = std::min(min_delta, lb - pruned_k);
                deltas.push_back(lb - pruned_k);
            };

            auto call_cb = [&](Cost current_k) {
                ++num_calls;
                if (current_k < k_min) {
                    deltas.push_back(k_min - current_k);
                    return true;
                }
                return false;
            };


            auto start = std::chrono::steady_clock::now();


            solved = editor.edit(k, result_cb, prune_cb, call_cb);


            auto end = std::chrono::steady_clock::now();
            auto duration = end - start;

            stats.emplace_back(k, editor.stats().allCalls(), std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count());

            if (std::chrono::seconds(0) <= timelimit && duration > timelimit)
                return Result::Unsolved();


            if (solved)
                continue;


            assert(!deltas.empty());
            size_t index = std::clamp<size_t>(quantile * deltas.size(), 0, deltas.size() - 1);
            assert(index < deltas.size());
            std::nth_element(deltas.begin(), deltas.begin() + static_cast<long>(index), deltas.end(), std::less<>());
            Cost quantile_delta = std::max(deltas[index], 1);
            assert(quantile_delta > 0);

            if (config.verbosity) {
                std::cout << "edit(" << std::setw(10) << k << "):";
                std::cout << " min_delta = " << std::setw(11) << min_delta;
                std::cout << " deltas size = " << std::setw(10) << deltas.size();
                std::cout << " quantile_delta = " << std::setw(10) << quantile_delta;
                std::cout << " num_calls = " << num_calls;
                std::cout << "\n";
            }

            k += quantile_delta;
        } while (!solved);

        return Result::Solved(solutions);
    }

    static Result search_incremental(const Instance &instance, const Configuration &config, std::vector<Stat> &stats, Cost delta) {
        if (delta < 1) {
            std::cerr << "delta must be at least 1" << std::endl;
            abort();
        }
        Editor editor(instance.graph.copy(), instance.costs, config);

        std::chrono::seconds timelimit(config.timelimit);

        Cost k = editor.initial_lower_bound();
        Cost min_remaining_cost = 0;
        bool solved;
        std::vector<Solution> solutions;

        do {
            size_t num_calls = 0;
            auto result_cb = [&](const std::vector<VertexPair> &edits) {
                solutions.emplace_back(instance, edits);
                min_remaining_cost = k - solutions.back().cost;
            };

            auto prune_cb = [](Cost /*pruned_k*/, Cost /*lb*/) {};

            auto call_cb = [&](Cost current_k) {
                ++num_calls;
                return current_k < min_remaining_cost;
            };


            auto start = std::chrono::steady_clock::now();


            solved = editor.edit(k, result_cb, prune_cb, call_cb);


            auto end = std::chrono::steady_clock::now();
            auto duration = end - start;

            stats.emplace_back(k, editor.stats().allCalls(), std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count());
            if (std::chrono::seconds(0) <= timelimit && duration > timelimit)
                return Result::Unsolved();


            if (config.verbosity) {
                std::cout << "edit(" << std::setw(10) << k << "):";
                std::cout << " delta = " << std::setw(6) << delta;
                std::cout << " num_calls = " << num_calls << "\n";
            }


            k += delta;
        } while (!solved);

        return Result::Solved(solutions);
    }

    /**
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
     * @param instance
     * @param config
     * @return
     */
    static Result search_exponential(const Instance &instance, const Configuration &config, std::vector<Stat> &stats,
            size_t max_history_length = 3, double desired_calls_factor = 2, double max_step_factor = 4,
            double quantile = 0.5) {
        Editor editor(instance.graph.copy(), instance.costs, config);

        std::chrono::seconds timelimit(config.timelimit);

        Cost k_init = editor.initial_lower_bound();
        Cost delta = 0;
        Cost min_remaining_cost = 0; // When a solutions is found, prune branches early when they are already worse than the found solution.
        bool solved;
        std::vector<Solution> solutions;
        std::deque<Cost> ks;
        std::deque<size_t> num_calls;

        do {
            Cost k = k_init + delta;

            Cost min_delta = std::numeric_limits<Cost>::max();
            std::vector<Cost> deltas;
            ks.push_back(k);
            num_calls.push_back(0);

            auto result_cb = [&](const std::vector<VertexPair> &edits) {
                solutions.emplace_back(instance, edits);
                min_remaining_cost = k - solutions.back().cost;
            };

            auto prune_cb = [&](Cost pruned_k, Cost lb) {
                min_delta = std::min(min_delta, lb - pruned_k);
                deltas.push_back(lb - pruned_k);
            };

            auto call_cb = [&](Cost current_k) {
                ++num_calls.back();
                if (current_k < min_remaining_cost) {
                    deltas.push_back(min_remaining_cost - current_k);
                    return true;
                }
                return false;
            };

            auto start = std::chrono::steady_clock::now();


            solved = editor.edit(k, result_cb, prune_cb, call_cb);


            auto end = std::chrono::steady_clock::now();
            auto duration = end - start;

            stats.emplace_back(k, editor.stats().allCalls(), std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count());

            if (std::chrono::seconds(0) <= timelimit && duration > timelimit)
                return Result::Unsolved();

            if (solved)
                break;


            if (ks.size() > max_history_length) {
                ks.pop_front();
                num_calls.pop_front();
            }


            assert(!deltas.empty());
            size_t index = std::clamp<size_t>(quantile * deltas.size(), 0, deltas.size() - 1);
            assert(index < deltas.size());
            std::nth_element(deltas.begin(), deltas.begin() + static_cast<long>(index), deltas.end(), std::less<>());
            Cost quantile_delta = std::max(deltas[index], 1);
            assert(quantile_delta > 0);

            Cost prev_delta = delta;
            Cost lo = std::max(1, prev_delta + quantile_delta);
            Cost hi = std::max<Cost>(lo, std::floor(max_step_factor * delta));
            Cost next_delta = find_next_k(ks, num_calls, desired_calls_factor) - k_init;
            delta = std::clamp(next_delta, lo, hi);

            if (config.verbosity) {
                std::cout << "edit(" << std::setw(10) << k << "):";
                std::cout << " lo = " << std::setw(6) << lo;
                std::cout << " next_delta = " << std::setw(6) << next_delta;
                std::cout << " hi = " << std::setw(6) << hi;
                std::cout << " delta = " << std::setw(6) << delta;
                std::cout << " num_calls = " << num_calls.back() << "\n";
            }

        } while (!solved);

        return Result::Solved(solutions);
    }

    /**
     *
     * @param ks
     * @param calls
     * @param factor
     * @return
     */
    static Cost find_next_k(const std::deque<Cost> &ks, const std::deque<size_t> &calls, double factor) {
        std::vector<double> xs;
        for (auto k : ks)
            xs.push_back(k);

        std::vector<double> ys;
        for (auto call : calls)
            ys.push_back(std::log(static_cast<double>(call)));

        double x_mean = 0;
        for (auto x : xs)
            x_mean += x;
        x_mean /= xs.size();

        double y_mean = 0;
        for (auto y : ys)
            y_mean += y;
        y_mean /= ys.size();


        double f1 = 0, f2 = 0;
        for (size_t i = 0; i < xs.size(); ++i) {
            f1 += (xs[i] - x_mean) * (ys[i] - y_mean);
            f2 += (xs[i] - x_mean) * (xs[i] - x_mean);
        }
        double beta = f1 / f2;
        double alpha = y_mean - beta * x_mean;

        double next_y = std::log(factor * calls.back());
        double next_x = (next_y - alpha) / beta;

        return std::ceil(next_x);
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_FPTSOLVER_H
