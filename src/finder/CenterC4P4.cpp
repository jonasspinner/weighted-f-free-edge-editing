//
// Created by jonas on 03.07.19.
//


#include "CenterC4P4.h"

#include "FinderI.h"


namespace Finder {
    /**
     * Calls callback for all P_4's and C_4's.
     *
     * @param graph
     * @param callback
     * @return
     */
    bool CenterC4P4::find(const Graph &graph, SubgraphCallback callback) {
        return find(graph, callback, neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
    }

    /**
     * Calls callback for all P_4's and C_4's. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
     *
     * @param graph
     * @param forbidden
     * @param callback
     * @return
     */
    bool CenterC4P4::find(const Graph &graph, const Graph &forbidden, SubgraphCallback callback) {
        return find(graph, callback, neighbors(graph, forbidden), non_neighbors(graph, forbidden),
                    valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    /**
     * Calls callback for all P_4's and C_4's having both u and v as vertices.
     *
     * @param uv
     * @param graph
     * @param callback
     * @return
     */
    bool CenterC4P4::find_near(VertexPair uv, const Graph &graph, SubgraphCallback callback) {
        return find_near(uv, callback, neighbors(graph), non_neighbors(graph), valid_edge(graph),
                         valid_non_edge(graph));
    }

    /**
     * Calls callback for all P_4's and C_4's having both u and v as vertices. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
     *
     * @param uv
     * @param graph
     * @param forbidden
     * @param callback
     * @return
     */
    bool CenterC4P4::find_near(VertexPair uv, const Graph &graph, const Graph &forbidden, SubgraphCallback callback) {
        return find_near(uv, callback, neighbors(graph, forbidden), non_neighbors(graph, forbidden),
                         valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    bool CenterC4P4::find_with_duplicates(const Graph &graph, const SubgraphCallback &callback) {
        return find_with_duplicates(graph, callback, neighbors(graph), non_neighbors(graph),
                                    valid_edge(graph), valid_non_edge(graph));
    }

    bool
    CenterC4P4::find_with_duplicates(const Graph &graph, const Graph &forbidden, const SubgraphCallback &callback) {
        return find_with_duplicates(graph, callback, neighbors(graph, forbidden), non_neighbors(graph, forbidden),
                                    valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    bool
    CenterC4P4::find_near_with_duplicates(VertexPair uv, const Graph &graph, const Graph &forbidden,
                                          const SubgraphCallback &callback) {
        // /*
        // version 1

        auto x = [&](Vertex a, Vertex b) -> size_t {
            return forbidden.hasEdge({a, b}) ? 1 : 0;
        };

        auto ensure_orientation = [](Subgraph &&subgraph) -> Subgraph {
            if (subgraph[0] > subgraph[3]) {
                std::swap(subgraph[0], subgraph[3]);
                std::swap(subgraph[1], subgraph[2]);
            }
            return std::move(subgraph);
        };

        return find_near(uv, graph, [&](Subgraph &&subgraph) {
            auto a = subgraph[0], b = subgraph[1], c = subgraph[2], d = subgraph[3];
            if (graph.hasEdge({a, d})) {
                // C4

                if (x(a, c) || x(b, d))
                    return false;

                if (x(a, b) + x(b, c) + x(c, d) + x(d, a) == 0) { // ab, bc, cd, ad
                    auto s1 = ensure_orientation(Subgraph{a, b, c, d}); // ad
                    assert(s1[0] < s1[3]);
                    if (callback(std::move(s1)))
                        return true;

                    auto s2 = ensure_orientation(Subgraph{b, c, d, a}); // ab
                    assert(s2[0] < s2[3]);
                    if (callback(std::move(s2)))
                        return true;

                    auto s3 = ensure_orientation(Subgraph{c, d, a, b}); // bc
                    assert(s3[0] < s3[3]);
                    if (callback(std::move(s3)))
                        return true;

                    auto s4 = ensure_orientation(Subgraph{d, a, b, c}); // cd
                    assert(s4[0] < s4[3]);
                    if (callback(std::move(s4)))
                        return true;

                } else if (x(a, b) && x(b, c) + x(c, d) + x(a, d) == 0) { // ab
                    auto s = ensure_orientation(Subgraph{b, c, d, a});
                    assert(s[0] < s[3]);
                    if (callback(std::move(s)))
                        return true;
                } else if (x(b, c) && x(a, b) + x(c, d) + x(a, d) == 0) { // bc
                    auto s = ensure_orientation(Subgraph{c, d, a, b});
                    assert(s[0] < s[3]);
                    if (callback(std::move(s)))
                        return true;
                } else if (x(c, d) && x(a, b) + x(b, c) + x(a, d) == 0) { // cd
                    auto s = ensure_orientation(Subgraph{d, a, b, c});
                    assert(s[0] < s[3]);
                    if (callback(std::move(s)))
                        return true;
                } else if (x(a, d) && x(a, b) + x(b, c) + x(c, d) == 0) { // ad
                    auto s = ensure_orientation(Subgraph{a, b, c, d});
                    assert(s[0] < s[3]);
                    if (callback(std::move(s)))
                        return true;
                }
            } else {
                // P4
                if (x(a, b) + x(a, c) + x(b, c) + x(b, d) + x(c, d) == 0) {
                    auto s = ensure_orientation(std::move(subgraph));
                    assert(s[0] < s[3]);
                    if (callback(std::move(s)))
                        return true;
                }
            }
            return false;
        });
        //*/


        /*
        // version 2
        const auto _neighbors = neighbors(graph);
        const auto _non_neighbors = non_neighbors(graph, forbidden);
        auto covered = [&](Vertex a, Vertex b) -> size_t {
            return forbidden.hasEdge({a, b}) ? 1 : 0;
        };


        bool conversion_edit_covered = forbidden.hasEdge(uv);

        if (graph.hasEdge(uv)) {
            // C4
            // u-v-a-b-u
            // b-u-v-a-b
            // a-b-u-v-a
            // v-a-b-u-v

            // P4
            // u-v-a-b
            // b-u-v-a
            // a-b-u-v


            // u-v-a
            A = _neighbors(v) & _non_neighbors(u);
            for (auto a : Graph::iterate(A)) {
                assert(!forbidden.hasEdge({u, a}));  // non_neighbors(graph, forbidden)
                if (forbidden.hasEdge({u, a}))
                    continue;
                if (forbidden.hasEdge({v, a})) {
                    if (conversion_edit_covered)
                        continue;
                    assert(!conversion_edit_covered);
                    conversion_edit_covered = true;
                }

                // u-v-a-b
                // u-v-a-b-u
                B = _neighbors(a) & _non_neighbors(v);
                for (auto b : Graph::iterate(B)) {
                    assert(!forbidden.hasEdge({v, b}));  // non_neighbors(graph, forbidden)
                    if (forbidden.hasEdge({v, b}))
                        continue;
                    if (forbidden.hasEdge({a, b})) {
                        if (conversion_edit_covered)
                            continue;
                        assert(!conversion_edit_covered);
                        conversion_edit_covered = true;
                    }

                    if (!conversion_edit_covered || !forbidden.hasEdge({u, b})) {
                        assert(graph.hasEdge({u, v})); assert(!graph.hasEdge({u, a})); assert(graph.hasEdge({v, a})); assert(!graph.hasEdge({v, b})); assert(graph.hasEdge({a, b}));
                        size_t _n_cov = 0; _n_cov += forbidden.hasEdge({u, v}); _n_cov += forbidden.hasEdge({u, a}); _n_cov += forbidden.hasEdge({u, b}); _n_cov += forbidden.hasEdge({v, a}); _n_cov += forbidden.hasEdge({v, b}); _n_cov += forbidden.hasEdge({a, b});
                        assert(_n_cov <= 1);
                        if (callback(Subgraph{u, v, a, b})) return true;
                    }


                    if (forbidden.hasEdge({a, b})) {
                        assert(conversion_edit_covered);
                        conversion_edit_covered = false;
                        assert(!forbidden.hasEdge(uv) && !forbidden.hasEdge({v, b}));
                    }
                }


                if (forbidden.hasEdge({v, a})) {
                    assert(conversion_edit_covered);
                    conversion_edit_covered = false;
                    assert(!forbidden.hasEdge(uv));
                }
            }


            // b-u-v-a
            // b-u-v-a-b
            B = _neighbors(u) & _non_neighbors(v);
            for (auto a : Graph::iterate(A)) {
                assert(!forbidden.hasEdge({u, a}));
                if (forbidden.hasEdge({v, a})) {
                    if (conversion_edit_covered)
                        continue;
                    assert(!conversion_edit_covered);
                    conversion_edit_covered = true;
                }

                for (auto b : Graph::iterate(B)) {
                    assert(!forbidden.hasEdge({v, b}));
                    if (forbidden.hasEdge({a, b})) {
                        if (conversion_edit_covered)
                            continue;
                        assert(!conversion_edit_covered);
                        conversion_edit_covered = true;
                    }


                    if (!conversion_edit_covered || !forbidden.hasEdge({a, b})) {
                        assert(graph.hasEdge({b, u})); assert(!graph.hasEdge({b, v})); assert(graph.hasEdge({u, v})); assert(!graph.hasEdge({u, a})); assert(graph.hasEdge({v, a}));
                        size_t _n_cov = 0; _n_cov += forbidden.hasEdge({u, v}); _n_cov += forbidden.hasEdge({u, a}); _n_cov += forbidden.hasEdge({u, b}); _n_cov += forbidden.hasEdge({v, a}); _n_cov += forbidden.hasEdge({v, b}); _n_cov += forbidden.hasEdge({a, b});
                        assert(_n_cov <= 1);
                        if (callback(Subgraph{u, v, a, b})) return true;
                    }

                    if (forbidden.hasEdge({a, b})) {
                        assert(conversion_edit_covered);
                        conversion_edit_covered = false;
                        assert(!forbidden.hasEdge(uv) && !forbidden.hasEdge({v, b}));
                    }
                }

                if (forbidden.hasEdge({v, a})) {
                    assert(conversion_edit_covered);
                    conversion_edit_covered = false;
                    assert(!forbidden.hasEdge(uv));
                }
            }


            // a-b-u-v
            // a-b-u-v-a
            B = _neighbors(u) & _non_neighbors(v);
            for (auto b : Graph::iterate(B)) {


                A = _neighbors(b) & _non_neighbors(u);
                for (auto a : Graph::iterate(A)) {

                }
            }
        } else {
            // C4
            // u-a-v-b-u
            // b-u-a-v-b

            // P4
            // u-a-v-b
            // b-u-a-v
            // u-a-b-v


            // u-a-v
            A = _neighbors(u) & _neighbors(v);
            for (auto a : Graph::iterate(A)) {
                // u-a-v-b
                // u-a-v-b-u
                B = _neighbors(v) & _non_neighbors(a);

                // b-u-a-v
                // b-u-a-v-b
                B = _neighbors(u) & _non_neighbors(a);
            }

            // u v
            A = _neighbors(u) & _non_neighbors(v);
            B = _neighbors(v) & _non_neighbors(u);
            for (auto a : Graph::iterate(A)) {
                for (auto b : Graph::iterate(B)) {
                    if (graph.hasEdge({a, b})) {
                        // u-a-b-v
                    }
                }
            }
        }
        */

        /*
        // version 3
        if (graph.hasEdge(uv)) {
            // u-v
            for (auto[a, b] : graph.edges()) {
                if (a == u || a == v || b == u || b == v)
                    continue;
                if (graph.hasEdge({u, a})) {

                } else if (graph.hasEdge({v, a})) {

                }
            }
        } else {
            // u-a b-v
            A = _neighbors(u) & _non_neighbors(v);
            B = _neighbors(v) & _non_neighbors(u);
            for (auto a : Graph::iterate(A)) {
                for (auto b : Graph::iterate(B)) {
                    if (graph.hasEdge({a, b})) {
                        // u-a-b-v
                    }
                }
            }
        }
        */

        return false;
    }

    bool CenterC4P4::for_all_conversionless_edits(const Subgraph &subgraph, const VertexPairCallback &callback) const {
        assert(subgraph.size() == 4);
        Vertex a = subgraph[0], b = subgraph[1], c = subgraph[2], d = subgraph[3];

        // Case P_4: +{a, d} results in a C_4
        // Case C_4: -{a, d} results in a P_4

        if (callback({a, b})) return true;
        if (callback({a, c})) return true;
        if (callback({b, c})) return true;
        if (callback({b, d})) return true;
        return callback({c, d});
    }

    void CenterC4P4::to_yaml(YAML::Emitter &out) const {
        using namespace YAML;
        out << BeginMap;
        out << Key << "name" << Value << "CenterC4P4";
        out << Key << "forbidden_subgraphs" << Value << Options::FSG::C4P4;
        out << EndMap;
    }

    template<typename F, typename G, typename H, typename I>
    bool
    CenterC4P4::find(const Graph &graph, const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge,
                     I valid_non_edge) {

        for (Vertex u : graph.vertices()) {
            V = neighbors(u);
            for (Vertex v : Graph::iterate(V)) {

                A = neighbors(u) & non_neighbors(v);
                B = neighbors(v) & non_neighbors(u);

                for (Vertex a : Graph::iterate(A)) {
                    for (Vertex b : Graph::iterate(B)) {

                        // assert that {a, u, v, b} is either a P_4 or a C_4
                        assert(valid_edge({a, u}));
                        assert(valid_non_edge({a, v})); /*assert(valid({a, b}));*/ assert(valid_edge({u, v}));
                        assert(valid_non_edge({u, b}));
                        assert(valid_edge({v, b}));

                        if (valid_non_edge({a, b}) && a < b) {
                            // P_4
                            if (callback(Subgraph{a, u, v, b})) return true;
                        }
                        if (valid_edge({a, b}) && a < std::min({u, v, b}) && u < b) {
                            // C_4
                            if (callback(Subgraph{a, u, v, b})) return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    template<typename F, typename G, typename H, typename I>
    bool
    CenterC4P4::find_near(VertexPair uv, const SubgraphCallback &callback, F neighbors, G non_neighbors, H valid_edge,
                          I valid_non_edge) {

        auto[u, v] = uv;
        assert(u < v);

        if (valid_non_edge(uv)) {
            // case P_4 u, v are at positions 1 and 4 in the path
            A = neighbors(u) & non_neighbors(v);
            B = neighbors(v) & non_neighbors(u);

            for (Vertex a : Graph::iterate(A)) {
                for (Vertex b : Graph::iterate(B)) {
                    if (valid_edge({a, b})) {

                        // assert that {u, a, b, v} is a P_4
                        assert(valid_edge({u, a}));
                        assert(valid_non_edge({u, b}));
                        assert(valid_non_edge({u, v}));
                        assert(valid_edge({a, b}));
                        assert(valid_non_edge({a, v}));
                        assert(valid_edge({b, v}));

                        // P_4
                        if (callback(Subgraph{u, a, b, v})) return true;
                    }
                }
            }

            // case P_4, C_4 u, v are at positions 1 and 3 in the path
            A = neighbors(u) & neighbors(v);
            for (Vertex a : Graph::iterate(A)) {
                B = neighbors(v) & non_neighbors(a);
                for (Vertex b : Graph::iterate(B)) {

                    // assert that {u, a, v, b} is either a C_4 or P_4
                    assert(valid_edge({u, a}));
                    assert(valid_non_edge({u, v})); /*assert(valid({u, b}));*/ assert(valid_edge({a, v}));
                    assert(valid_non_edge({a, b}));
                    assert(valid_edge({v, b}));

                    if (valid_edge({u, b}) && u < std::min({a, v, b}) && a < b) {
                        // C_4
                        if (callback(Subgraph{u, a, v, b})) return true;
                    } else if (valid_non_edge({u, b})) {
                        // P_4
                        if (callback(Subgraph{u, a, v, b})) return true;
                    }
                }
            }

            // case P_4, C_4 u, v are at positions 2 and 4 in the path
            B = neighbors(u) & neighbors(v);
            for (Vertex b : Graph::iterate(B)) {
                A = neighbors(u) & non_neighbors(b);
                for (Vertex a : Graph::iterate(A)) {

                    // assert that {a, u, b, b} is either a C_4 or P_4
                    assert(valid_edge({a, u}));
                    assert(valid_non_edge({a, b})); /*assert(valid({a, v}));*/ assert(valid_edge({u, b}));
                    assert(valid_non_edge({u, v}));
                    assert(valid_edge({b, v}));

                    if (valid_edge({a, v}) && a < std::min({u, b, v})) {
                        // C_4
                        if (callback(Subgraph{a, u, b, v})) return true;
                    } else if (valid_non_edge({a, v})) {
                        // P_4
                        if (callback(Subgraph{a, u, b, v})) return true;
                    }
                }
            }

        } else if (valid_edge(uv)) {
            // case C_4 u, v are at positions 1 and 4 in the path
            A = neighbors(u) & non_neighbors(v);
            B = neighbors(v) & non_neighbors(u);

            for (Vertex a : Graph::iterate(A)) {
                for (Vertex b : Graph::iterate(B)) {
                    if (valid_edge({a, b}) && u < std::min({a, b, v}) && a < v) {

                        // assert that {u, a, b, v} is a C_4
                        assert(valid_edge({u, a}));
                        assert(valid_non_edge({u, b}));
                        assert(valid_edge({u, v}));
                        assert(valid_edge({a, b}));
                        assert(valid_non_edge({a, v}));
                        assert(valid_edge({b, v}));

                        // C_4
                        if (callback(Subgraph{u, a, b, v})) return true;
                    }
                }
            }

            // case C_4, P_4 u, v are at positions 2 and 3 in the path
            A = neighbors(u) & non_neighbors(v);
            B = neighbors(v) & non_neighbors(u);

            for (Vertex a : Graph::iterate(A)) {
                for (Vertex b : Graph::iterate(B)) {

                    // assert that {a, u, v, b} is either a C_4 or a P_4
                    assert(valid_edge({a, u}));
                    assert(valid_non_edge({a, v})); /*assert(valid({a, b}));*/ assert(valid_edge({u, v}));
                    assert(valid_non_edge({u, b}));
                    assert(valid_edge({v, b}));

                    if (valid_edge({a, b}) && a < std::min({u, v, b})) {
                        // C_4
                        if (callback(Subgraph{a, u, v, b})) return true;
                    } else if (valid_non_edge({a, b})) {
                        // P_4
                        if (callback(Subgraph{a, u, v, b})) return true;
                    }
                }
            }

            // case C_4, P_4 u, v are at positions 1 and 2 in the path
            A = neighbors(v) & non_neighbors(u);

            for (Vertex a : Graph::iterate(A)) {
                B = neighbors(a) & non_neighbors(v);
                B[u] = false;

                for (Vertex b : Graph::iterate(B)) {

                    // assert that {u, v, a, b} is either a C_4 or a P_4
                    assert(valid_edge({u, v}));
                    assert(valid_non_edge({u, a})); /*assert(valid({u, b}));*/ assert(valid_edge({v, a}));
                    assert(valid_non_edge({v, b}));
                    assert(valid_edge({a, b}));

                    if (valid_edge({u, b}) && u < std::min({v, a, b}) && v < b) {
                        // C_4
                        if (callback(Subgraph{u, v, a, b})) return true;
                    } else if (valid_non_edge({u, b})) {
                        // P_4
                        if (callback(Subgraph{u, v, a, b})) return true;
                    }
                }
            }

            // case C_4, P_4 u, v are at positions 3 and 4 in the path
            A = neighbors(u) & non_neighbors(v);

            for (Vertex a : Graph::iterate(A)) {
                B = neighbors(a) & non_neighbors(u);
                B[v] = false;

                for (Vertex b : Graph::iterate(B)) {

                    // assert that {b, a, u, v} is either a C_4 or P_4
                    assert(valid_edge({b, a}));
                    assert(valid_non_edge({b, u})); /*assert(valid({b, v}));*/ assert(valid_edge({a, u}));
                    assert(valid_non_edge({a, v}));
                    assert(valid_edge({u, v}));

                    if (valid_edge({b, v}) && b < std::min({a, u, v})) {
                        // C_4
                        if (callback(Subgraph{b, a, u, v})) return true;
                    } else if (valid_non_edge({b, v})) {
                        // P_4
                        if (callback(Subgraph{b, a, u, v})) return true;
                    }
                }
            }
        }
        return false;
    }

    template<typename F, typename G, typename H, typename I>
    bool CenterC4P4::find_with_duplicates(const Graph &graph, const FinderI::SubgraphCallback &callback, F neighbors,
                                          G non_neighbors, H valid_edge, I valid_non_edge) {

        auto ensure_orientation = [](Subgraph &&subgraph) -> Subgraph {
            if (subgraph[0] > subgraph[3]) {
                std::swap(subgraph[0], subgraph[3]);
                std::swap(subgraph[1], subgraph[2]);
            }
            return std::move(subgraph);
        };

        for (auto uv : graph.edges()) {
            if (!valid_edge(uv))
                continue;
            auto[u, v] = uv;

            A = neighbors(u) & non_neighbors(v);
            B = neighbors(v) & non_neighbors(u);

            for (auto a : Graph::iterate(A)) {
                for (auto b : Graph::iterate(B)) {
                    // assert that {a, u, v, b} is either a P_4 or a C_4
                    assert(valid_edge({a, u}));
                    assert(valid_non_edge({a, v})); /*assert(valid({a, b}));*/ assert(valid_edge({u, v}));
                    assert(valid_non_edge({u, b}));
                    assert(valid_edge({v, b}));
                    // C_4 or P_4
                    auto s = ensure_orientation(Subgraph{a, u, v, b});
                    assert(s[0] < s[3]);
                    if (callback(std::move(s))) return true;
                }
            }
        }
        return false;
    }
}
