//
// Created by jonas on 27.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_PRUNEDDELTA_H
#define WEIGHTED_F_FREE_EDGE_EDITING_PRUNEDDELTA_H


/**
 * From initial implementation:
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

size_t index = std::clamp<size_t>(quantile * deltas.size(), 0, deltas.size() - 1);
std::nth_element(deltas.begin(), deltas.begin() + static_cast<long>(index), deltas.end(), std::less<>());
Cost quantile_delta = std::max(deltas[index], 1);

k += quantile_delta;
 */


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
 */
class PrunedDelta : public SearchStrategyI {
    Cost m_current_search_k{-1};
    std::vector<Cost> m_prune_deltas;
    float m_quantile;
public:

    void set_initial_search_k(Cost initial_search_k) override {
        m_current_search_k = initial_search_k;
    }

    std::optional<Cost> get_next_search_k() override {
        Cost delta = get_quantile(m_prune_deltas, m_quantile);

        m_current_search_k += delta;

        return m_current_search_k;
    }

    void register_call(Cost /*k*/) override {}

    void register_bound(Cost k, Cost lower_bound) override {
        if (lower_bound == invalid_cost) return;

        m_prune_deltas.push_back(lower_bound - k);
    }

private:
    static Cost get_quantile(std::vector<Cost> &values, float quantile) {
        assert(!values.empty());

        size_t index = std::clamp<size_t>(quantile * values.size(), 0, values.size() - 1);
        assert(index < values.size());

        std::nth_element(values.begin(), values.begin() + static_cast<long>(index), values.end(), std::less<>());

        return values[index];
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_PRUNEDDELTA_H
