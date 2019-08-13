//
// Created by jonas on 12.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LB_ARW_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LB_ARW_H


#include <utility>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <iostream>
#include <memory>
#include <random>

#include "../config.hpp"

#include "../Options.hpp"
#include "../Finder/Finder.hpp"
#include "../Finder/SubgraphStats.hpp"
#include "../LowerBound/Lower_Bound.hpp"
#include "../Graph/ValueMatrix.hpp"
#include "../util.hpp"


#include "../interfaces/LowerBoundI.h"
#include "../graph/VertexPairMap.h"
#include "../consumer/SubgraphStats.h"

namespace Consumer
{
    class ARW : LowerBoundI
    {
        
    public:
        static constexpr char const *name = "ARW";
        class Lower_Bound_Storage_type {
        private:
            std::vector<Subgraph> bound;
        public:
            template <typename iterator_t>
            void add(const iterator_t& begin, const iterator_t& end)
            {
                bound.emplace_back(begin, end);
            }

            void add(const Subgraph& sg)
            {
                bound.emplace_back(sg);
            }

            /**
             * Removes any forbidden subgraphs that contain the given vertex pair (u, v).
             */
            void remove(const Graph& graph, const VertexPairMap<bool>& edited, Vertex u, Vertex v)
            {
                if (edited[{u, v}])
                {
                    abort();
                }
                auto new_end = std::remove_if(bound.begin(), bound.end(), [&](const Subgraph& subgraph) -> bool {
                    auto vertexPairs = subgraph.vertexPairs();
                    return std::any_of(vertexPairs.begin(), vertexPairs.end(), [&](VertexPair xy) {
                        if (!edited[xy]) {
                            return xy == VertexPair{u, v};
                        }
                        return false;
                    });
                });

                bound.erase(new_end, bound.end());
            }

            void clear()
            {
                bound.clear();
            }


            [[nodiscard]] const std::vector<Subgraph>& get_bound() const
            {
                return bound;
            }

            std::vector<Subgraph>& get_bound()
            {
                return bound;
            }

            const Subgraph& operator[](size_t i) const
            {
                return bound[i];
            }

            [[nodiscard]] std::vector<Vertex> as_vector(size_t i) const
            {
                return std::vector<Vertex>(bound[i].begin(), bound[i].end());
            }

            [[nodiscard]] size_t size() const
            {
                return bound.size();
            }

            [[nodiscard]] bool empty() const
            {
                return bound.empty();
            }


            void assert_valid(const Graph&
#ifndef NDEBUG
            graph
#endif
                    , const VertexPairMap<bool>&
#ifndef NDEBUG
            edited
#endif
            ) const
            {
#ifndef NDEBUG
                VertexPairMap<bool> in_bound(graph.size());

                for (const auto& fs : bound)
                {
                    for (VertexPair uv : fs.vertexPairs()) {
                        if (!edited[uv]) {
                            assert(!in_bound[uv]);
                            in_bound[uv];
                        }
                    }
                }
#endif
            }

        };
        // using Lower_Bound_Storage_type = ::Lower_Bound::Lower_Bound<Mode, Restriction, Conversion, Graph, VertexPairMap<bool>, length>;
        //using Subgraph_Stats_type = ::Finder::Subgraph_Stats<Finder_impl, Graph, VertexPairMap<bool>, Mode, Restriction, Conversion, length>;
        //using Subgraph = typename Lower_Bound_Storage_type::Subgraph;

        static constexpr bool needs_subgraph_stats = true;

        struct State : public StateI {
            Lower_Bound_Storage_type lb;
            std::unique_ptr<StateI> copy() override {
                return std::make_unique<State>(*this);
            }
        };
    private:

        const Graph &m_graph;
        const VertexPairMap<bool> &m_edited;
        VertexPairMap<bool> candidate_pairs_used;
        Graph bound_uses;
        SubgraphStats subgraph_stats;

    public:
        ARW(std::shared_ptr<FinderI> finder_ptr, const Instance &instance, const VertexPairMap<bool> &marked) : LowerBoundI(std::move(finder_ptr)), m_graph(instance.graph), m_edited(marked), candidate_pairs_used(instance.graph.size()), bound_uses(instance.graph.size()), subgraph_stats(finder_ptr, instance, marked) {}

        std::unique_ptr<StateI> initialize(Cost k) override
        {
            auto ptr = std::make_unique<State>();
            State& state = *ptr;

            bound_uses.clear_edges();

            finder->find(bound_uses, [&](const Subgraph& path)
            {
                bool touches_bound = false;
                for (VertexPair uv : path.vertexPairs()) {
                    if (!m_edited[uv] && bound_uses.has_edge(uv)) {
                        touches_bound = true;
                    }
                }

                if (!touches_bound)
                {
                    state.lb.add(path);

                    for (VertexPair uv : path.vertexPairs())
                        if (!m_edited[uv])
                            bound_uses.set_edge(uv);
                }

                // Assumption: if the bound is too high, initialize will be called again anyway.
                return state.lb.size() > k;
            });

            return ptr;
        }

        void before_mark_and_edit(StateI& state_, VertexPair uv) override
        {
            auto state = dynamic_cast<State &>(state_);

            std::vector<Subgraph>& lb = state.lb.get_bound();

            for (size_t i = 0; i < lb.size();)
            {
                bool has_uv = false;
                for (VertexPair xy : lb[i].vertexPairs()) {
                    if (!m_edited[xy]) {
                        if (uv == xy) has_uv = true;
                    }
                }

                if (has_uv)
                {
                    lb[i] = lb.back();
                    lb.pop_back();
                }
                else
                {
                    ++i;
                }
            }
        }

        void after_mark_and_edit(StateI& state_, VertexPair uv) override
        {
            auto state = dynamic_cast<State &>(state_);
            initialize_bound_uses(state, m_graph, m_edited);

            finder->find_near(uv, bound_uses, [&](const Subgraph& path)
            {
                bool touches_bound = false;
                for (VertexPair uv : path.vertexPairs()) {
                    if (!m_edited[uv]) {
                        if (bound_uses.has_edge(uv)) {
                            touches_bound = true;
                        }
                    }
                }

                if (!touches_bound)
                {
                    state.lb.add(path);

                    for (VertexPair uv : path.vertexPairs()) {
                        if (!m_edited[uv]) {
                            bound_uses.set_edge(uv);
                        }
                    }
                }

                return false;
            });
        }

        void before_mark(StateI&, VertexPair) override
        {
        }

        void after_mark(StateI& state_, VertexPair uv) override
        {
            after_mark_and_edit(state_, uv);
        }

        Cost result(StateI& state_, Cost k) override
        {
            auto state = dynamic_cast<State &>(state_);
            if (state.lb.size() <= k)
            {
                initialize_bound_uses(state, m_graph, m_edited);
                find_lb_2_improvements_v2(state, subgraph_stats, k, m_graph, m_edited);
            }

            return state.lb.size();
        }

    private:
        void initialize_bound_uses(const State& state, const Graph& graph, const VertexPairMap<bool> &edited)
        {
            bound_uses.clear_edges();

            for (const Subgraph& path : state.lb.get_bound())
            {
                for (VertexPair uv : path.vertexPairs()) {
                    if (!edited[uv]) {
                        bound_uses.set_edge(uv);
                    }
                }
            }
        }


        void find_lb_2_improvements_v2(State& state, const SubgraphStats& subgraph_stats, size_t k, const Graph &g, const VertexPairMap<bool> &e) {
            std::vector<Subgraph>& lb = state.lb.get_bound();

            std::mt19937_64 gen(42 * subgraph_stats.num_subgraphs + subgraph_stats.sum_subgraphs_per_edge);

            bool improvement_found;
            bool bound_changed;
            size_t rounds_no_improvement = 0;
            do
            {
                improvement_found = false;
                bound_changed = false;

                std::shuffle(lb.begin(), lb.end(), gen);

                for (size_t fsi = 0; fsi < lb.size(); ++fsi)
                {
                    bool early_exit = find_partner(state, subgraph_stats, k, g, e, fsi, bound_changed, improvement_found, gen);
                    if (early_exit) return;
                }

                if (improvement_found)
                {
                    rounds_no_improvement = 0;
                }
                else
                {
                    ++rounds_no_improvement;
                }
            } while (improvement_found || (rounds_no_improvement < 5 && bound_changed));

        }

        bool find_partner(State& state, const SubgraphStats& subgraph_stats, size_t k, const Graph &g, const VertexPairMap<bool> &e, size_t fsi, bool& bound_changed, bool& improvement_found, std::mt19937_64 &gen) {
            std::uniform_real_distribution<float> prob(.0, 1.0);

            std::vector<Subgraph>& lb = state.lb.get_bound();

            auto fs = lb[fsi];

            // See how many candidates we have
            // size_t num_pairs = 0, num_neighbors = 0;
            const auto [num_pairs, num_neighbors] = count_neighbors(subgraph_stats, g, e, fs);
            const auto pairs = get_pairs(g, e, fs);
            // const auto [num_pairs, num_neighbors, pairs] = get_pairs(subgraph_stats, g, e, fs);

            // Skip if this forbidden subgraph if it is the only one we can add.
            if (num_pairs > 1) {
                // Remove fs from lower bound
                remove_from_bound_uses(pairs, bound_uses);
                add_to_candidate_pairs_used(pairs, candidate_pairs_used);

                // Collect candidates
                const auto candidates_per_pair = get_candidates(finder, g, e, candidate_pairs_used, pairs, bound_uses);

                // Remove fs again from lower bound
                remove_from_bound_uses(pairs, bound_uses);
                remove_from_candidate_pairs_used(pairs, candidate_pairs_used);

                const bool random_switch = prob(gen) < 0.3;
                size_t min_candidate_neighbors = num_neighbors;
                size_t num_candidates_considered = 0;

                auto min_candidate = fs;
                size_t min_pairs = num_pairs;

                bool found_partner = false;

                for (size_t pi = 0; pi < pairs.size(); ++pi) {
                    for (auto cand_fs : candidates_per_pair[pi]) {
                        assert(!found_partner);

                        add_to_bound_uses(g, e, cand_fs, bound_uses);
                        const auto [cand_pairs, cand_neighbors] = count_neighbors(subgraph_stats, g, e, cand_fs);

                        ++num_candidates_considered;

                        if (cand_pairs == 1 || (min_pairs > 1 && (
                                (!random_switch && cand_neighbors < min_candidate_neighbors) ||
                                (random_switch && prob(gen) < 1.0 / num_candidates_considered)))) {
                            min_pairs = cand_pairs;
                            min_candidate = cand_fs;
                            min_candidate_neighbors = cand_neighbors;
                        }

                        // Look only for later pairs - combinations with subgraphs listed for earlier pairs have already been considered!
                        for (size_t ppi = pi + 1; ppi < pairs.size(); ++ppi) {
                            const auto partner_pair = pairs[ppi];

                            if (bound_uses.has_edge({partner_pair.first, partner_pair.second})) continue;

                            for (auto partner_fs : candidates_per_pair[ppi]) {
                                bool touches_bound = is_Subgraphouching_bound(g, e, partner_fs, bound_uses);

                                if (!touches_bound) {
                                    // Success!
                                    found_partner = true;
                                    improvement_found = true;

                                    // Directly add the partner to the lower bound, continue search to see if there is more than one partner
                                    lb.push_back(partner_fs);

                                    add_to_bound_uses(g, e, partner_fs, bound_uses);
                                }
                            }

                        }

                        if (found_partner) {
                            lb[fsi] = cand_fs;

                            state.lb.assert_valid(g, e);
                            break;
                        } else {
                            remove_from_bound_uses(g, e, cand_fs, bound_uses);
                        }
                    }

                    if (found_partner) break;
                }

                if (!found_partner) {
                    if (min_candidate != fs) {
                        add_to_bound_uses(g, e, min_candidate, bound_uses);

                        lb[fsi] = min_candidate;
                        state.lb.assert_valid(g, e);

                        bound_changed = true;
                    } else {
                        // Add fs back to lower bound
                        add_to_bound_uses(pairs, bound_uses);
                    }
                } else if (k < state.lb.size()) {
                    return true; //exit early
                }
            }

            return false; // normal exit
        }

        static std::vector<std::vector<Subgraph>> get_candidates(FinderI &finder, const Graph &g, const VertexPairMap<bool> &e, const VertexPairMap<bool> &candidate_pairs_used, const std::vector<std::pair<Vertex, Vertex>> &pairs, Graph &bound_uses) {

            std::vector<std::vector<Subgraph>> candidates_per_pair(pairs.size());

            for (size_t pi = 0; pi < pairs.size(); ++pi) {
                auto p = pairs[pi];

                auto cb = [&](const Subgraph &sg) {

                    size_t pairs_covered = 0;

                    for (VertexPair uv : sg.vertexPairs())
                        if (!e[uv])
                            pairs_covered += candidate_pairs_used[uv];

                    if (pairs_covered < pairs.size()) {
                        candidates_per_pair[pi].push_back(sg);
                    }

                    return false;
                };

                finder.find_near({p.first, p.second}, bound_uses, cb);


                // Set node pair to avoid getting the same candidates twice.
                // WARNING: this assumes all candidates are listed for all node pairs, i.e., clean_graph_structure has not been called!
                bound_uses.set_edge({p.first, p.second});
            }
            return candidates_per_pair;
        }

        static std::vector<std::pair<Vertex, Vertex>> get_pairs(const Graph &g, const VertexPairMap<bool> &e, const Subgraph& subgraph) {
            std::vector<std::pair<Vertex, Vertex>> pairs;
            for (VertexPair uv : subgraph.vertexPairs()) {
                if (!e[uv]) {
                    pairs.emplace_back(uv.u, uv.v);
                }
            }
            return pairs;
        }

        static std::tuple<size_t, size_t> count_neighbors(const SubgraphStats& subgraph_stats, const Graph &g, const VertexPairMap<bool> &e, const Subgraph& fs) {
            size_t num_pairs = 0, num_neighbors = 0;

            for (VertexPair uv : fs.vertexPairs()) {
                if (!e[uv]) {
                    size_t nn = subgraph_stats.num_subgraphs_per_edge[uv];
                    num_neighbors += nn;
                    if (nn > 1) ++num_pairs;
                }
            }
            return {num_pairs, num_neighbors};
        }

        static void add_to_bound_uses(const Graph &g, const VertexPairMap<bool> &e, const Subgraph& subgraph, Graph& bound_uses) {
            for (VertexPair uv : subgraph.vertexPairs()) {
                if (!e[uv]) {
                    assert(!bound_uses.has_edge(uv)); // check whether it holds
                    bound_uses.set_edge(uv);
                }
            }
        }

        static void remove_from_bound_uses(const Graph &g, const VertexPairMap<bool> &e, const Subgraph& subgraph, Graph& bound_uses) {
            for (VertexPair uv : subgraph.vertexPairs()) {
                if (!e[uv]) {
                    bound_uses.clear_edge(uv);
                }
            }
        }

        static bool is_subgraph_touching_bound(const Graph &g, const VertexPairMap<bool> &e, const Subgraph& subgraph, const Graph& bound_uses) {
            for (VertexPair uv : subgraph.vertexPairs()) {
                if (!e[uv]) {
                    if (bound_uses.has_edge(uv)) {
                        return true;
                    }
                }
            }
            return false;
        }


        static void remove_from_candidate_pairs_used(const std::vector<std::pair<Vertex, Vertex>> &pairs, VertexPairMap<bool> &candidate_pairs_used) {
            for (auto p : pairs) {
                candidate_pairs_used[{p.first, p.second}] = false;
            }
        }

        static void remove_from_bound_uses(const std::vector<std::pair<Vertex, Vertex>> &pairs, Graph &bound_uses) {
            for (auto p : pairs) {
                bound_uses.clear_edge({p.first, p.second});
            }
        }

        static void add_to_candidate_pairs_used(const std::vector<std::pair<Vertex, Vertex>> &pairs, VertexPairMap<bool> &candidate_pairs_used){
            for (auto p : pairs) {
                candidate_pairs_used[{p.first, p.second}] = true;
            }
        }

        static void add_to_bound_uses(const std::vector<std::pair<Vertex, Vertex>> &pairs, Graph& bound_uses) {
            for (auto p : pairs) {
                bound_uses.set_edge({p.first, p.second});
            }
        }

    };
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_LB_ARW_H
