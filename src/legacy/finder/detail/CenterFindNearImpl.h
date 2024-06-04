//
// Created by jonas on 31.07.19.
//


#ifndef WEIGHTED_F_FREE_EDGE_EDITING_CENTERFINDNEARIMPL_H
#define WEIGHTED_F_FREE_EDGE_EDITING_CENTERFINDNEARIMPL_H


#include "CenterFindImpl.h"


namespace detail {

    template <int k, bool with_cycles>
    class FindNearImpl {
        static_assert(k >= 4);
    public:
        /**
         * Find paths (and cycles) with length of at least 4 containing the vertices u and v. Each path is listed exactly once.
         *
         * The current implementation is inefficient. It lists all paths (and cycles) with length of at least 4 and filters them.
         *
         * @param graph
         * @param uv
         * @return
         */
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find_near(const Graph& graph, VertexPair uv, SubgraphCallback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            return CenterFindImpl<k, with_cycles>::find(graph, [&](Subgraph &&subgraph, Vertex) {
                auto vertices = subgraph.vertices();

                bool has_u = std::any_of(vertices.begin(), vertices.end(), [&](Vertex x) { return x == uv.u; });
                bool has_v = std::any_of(vertices.begin(), vertices.end(), [&](Vertex x) { return x == uv.v; });

                if (has_u && has_v) {
                    return callback(std::move(subgraph));
                } else {
                    return false;
                }
            }, neighbors, non_neighbors, valid_edge, valid_non_edge);
        }

        /**
         * Find paths (and cycles) with length of at least 4 containing the vertices u and v. Each path is listed exactly once.
         *
         * @param graph
         * @param uv
         * @return
         */
        /*
       template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
       static bool find_near(const Graph& graph, VertexPair uv, SubgraphCallback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
           auto [u, v] = uv;

           std::cout << "find_near " << uv << "\n";

           size_t d_init = 1, d_max = 0;
           if (valid_edge(uv)) {
               d_init = 0; d_max = 0;
           } else if (valid_non_edge(uv)) {
               d_init = 1; d_max = k - 2;
           }

           // For each possible distance between u and v on a path (or cycle).
           for (size_t d = d_init; d <= d_max; ++d) {
               std::cout << "d = " << d << "\n";
               std::vector<Vertex> init_path(d + 2);
               init_path[0] = u;
               init_path[d + 1] = v;

               // Fill path between u and v
               if (fill_inner(init_path, 0, d + 1, d + 2, graph, [&](const std::vector<Vertex> &vertices) {
                   assert(vertices.size() == d + 2);
                   std::vector<Vertex> path(k);

                   for (size_t i = 0; i < d + 2; ++i)
                       path[i] = vertices[i];

                   // Vertices is filled with a path. Move it on the result path
                   for (size_t offset = 0; offset < k - d - 1; ++offset) {
                       std::cout << "offset = " << offset << "\n";

                       for (size_t i = 0; i < d + 2; ++i)
                           path[i + offset] = vertices[i];

                       // Fill remaining part of the path
                       if (fill_outer(path, offset, offset + d + 1, graph, [&](Subgraph &&subgraph) {

                           if (valid_edge({subgraph[0], subgraph[k - 1]})) {
                               // C_k
#ifndef NDEBUG
                               assert(subgraph.size() == k);
                               for (unsigned i = 0; i < k; ++i)
                                   for (unsigned j = i + 1; j < k; ++j)
                                       if (j - i == 1 || j - i == k - 1) {
                                           assert(valid_edge({subgraph[i], subgraph[j]}));
                                       } else {
                                           assert(valid_non_edge({subgraph[i], subgraph[j]}));
                                       }
#endif
                               auto sg_vertices = subgraph.vertices();
                               Vertex min_vertex = *std::min_element(sg_vertices.begin(), sg_vertices.end());

                               if (subgraph[0] == min_vertex) {
                                   std::cout << "=> " << subgraph << "\n";
                                   if (callback(std::move(subgraph))) return true;
                               }
                           } else if (valid_non_edge({subgraph[0], subgraph[k - 1]})) {
                               // P_k
#ifndef NDEBUG
                               assert(subgraph.size() == k);
                               for (unsigned i = 0; i < k; ++i)
                                   for (unsigned j = i + 1; j < k; ++j)
                                       if (j - i == 1) {
                                           assert(valid_edge({subgraph[i], subgraph[j]}));
                                       } else {
                                           assert(valid_non_edge({subgraph[i], subgraph[j]}));
                                       }
#endif
                               std::cout << "=> " << subgraph << "\n";
                               if (callback(std::move(subgraph))) return true;
                           }
                           return false;
                       }, neighbors, non_neighbors, valid_edge, valid_non_edge)) return true;
                   }

                   // path [0..k-1]
                   // vertices [0..d+1]
                   for (size_t i = 0; i < d + 2; ++i)
                       path[i] = path[i + k - d - 2 ];

                   return false;
               }, neighbors, non_neighbors, valid_edge, valid_non_edge)) return true;
           }
           return false;
       }
       */

    private:
        template <typename Callback, typename F, typename G, typename H, typename I>
        static bool fill_inner(std::vector<Vertex> &vertices, size_t l, size_t r, size_t size, const Graph& graph, Callback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {

            std::cout << "fill_inner\t[ ";
            for (size_t i = 0; i <= l; ++i)
                std::cout << vertices[i] << " ";
            for (size_t i = l + 1; i < r; ++i)
                std::cout << "_ ";
            for (size_t i = r; i < size; ++i)
                std::cout << vertices[i] << " ";
            std::cout << "] l = " << l << ", r = " << r << "\n";


            if (r - l == 1) {
                if (valid_edge({vertices[l], vertices[r]}))
                    if (callback(vertices)) return true;
                return false;
            } else if (r - l == 2){
                // ajdacent to both p_l and p_r and not p_1, ..., p_{l-1}, p_{r+1}, ..., p_k
                auto C = neighbors(vertices[l]) & neighbors(vertices[r]);
                for (size_t i = 0;     i < l;    ++i) { C &= non_neighbors(vertices[i]); }
                for (size_t i = r + 1; i < size; ++i) { C &= non_neighbors(vertices[i]); }

                for (Vertex c : Graph::iterate(C)) {
                    vertices[l + 1] = c;
#ifndef NDEBUG
                    for (unsigned i = 0; i < size; ++i)
                        for (unsigned j = i + 1; j < size; ++j)
                            if (j - i == 1) {
                                assert(valid_edge({vertices[i], vertices[j]}));
                            } else if (j - i != k - 1){
                                assert(valid_non_edge({vertices[i], vertices[j]}));
                            }
#endif
                    if (callback(vertices)) return true;
                }
                return false;
            }
            /*
            auto C = graph.full_adjacency_row();
            for (size_t i = 0; i <= l; ++i)
                C[vertices[i]] = false;
            for (size_t i = r; i < size; ++i)
                C[vertices[i]] = false;
            for (size_t i = 0; i < l; ++i)
                C &= non_neighbors(vertices[i]);
            for (size_t i = r + 1; i < size; ++i)
                C &= non_neighbors(vertices[i]);

            auto A = C & neighbors(vertices[l]) & non_neighbors(vertices[r]);
            auto B = C & neighbors(vertices[r]) & non_neighbors(vertices[l]);
            */

            auto A = neighbors(vertices[l]); // adjacent to p_l, non adjacent to p_1, ..., p_{l-1}, p_r, ..., p_k
            for (size_t i = 0; i < l;    ++i) { A &= non_neighbors(vertices[i]); }
            for (size_t i = r; i < size; ++i) { A &= non_neighbors(vertices[i]); }

            auto B = neighbors(vertices[r]); // adjacent to p_r, non adjacent to p_1, ..., p_l, p_{r+1}, ..., p_k
            for (size_t i = 0;     i <= l;   ++i) { A &= non_neighbors(vertices[i]); }
            for (size_t i = r + 1; i < size; ++i) { A &= non_neighbors(vertices[i]); }

            for (Vertex a : Graph::iterate(A)) {
                for (Vertex b : Graph::iterate(B)) {
                    vertices[l + 1] = a;
                    vertices[r - 1] = b;

                    if (fill_inner(vertices, l + 1, r - 1, size, graph, callback, neighbors, non_neighbors, valid_edge, valid_non_edge)) return true;
                }
            }

            return false;
        }

        template <typename Callback, typename F, typename G, typename H, typename I>
        static bool fill_outer(std::vector<Vertex> &vertices, size_t l, size_t r, const Graph &graph, Callback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            assert(vertices.size() == k);

            std::cout << "fill_outer\t[ ";
            for (size_t i = 0; i < l; ++i)
                std::cout << "_ ";
            for (size_t i = l; i <= r; ++i)
                std::cout << vertices[i] << " ";
            for (size_t i = r + 1; i < vertices.size(); ++i)
                std::cout << "_ ";
            std::cout << "]\n";


            if (l == 0 && r == k - 1) {
                if (valid_edge({vertices.front(), vertices.back()}) || valid_edge({vertices.front(), vertices.back()})) {
                    Subgraph subgraph{};
                    for (Vertex u : vertices)
                        subgraph.push_back(u);

                    return callback(std::move(subgraph));
                }
            } else if (l == 0) {
                /*
                auto A = neighbors(vertices[r]);
                for (size_t i = l; i <= r; ++i)
                    A[vertices[i]] = false;

                if (r < vertices.size() - 2)
                    A &= non_neighbors(vertices[0]);

                for (size_t i = 1; i < r; ++i)
                    A &= non_neighbors(vertices[i]);
                */

                auto A = neighbors(vertices[r]); // adjacent to p_r, non adjacent to p_2, ..., p_{r-1}, if r != k - 1: non adjacent to p_1
                for (size_t i = 1; i < r; ++i)
                    A &= non_neighbors(vertices[i]);
                if (r != k - 2)
                    A &= non_neighbors(vertices[0]);

                for (Vertex a : Graph::iterate(A)) {
                    vertices[r + 1] = a;
                    if (fill_outer(vertices, 0, r + 1, graph, callback, neighbors, non_neighbors, valid_edge, valid_non_edge)) return true;
                }
            } else if (r == k - 1) {
                /*
                auto A = neighbors(vertices[l]);
                for (size_t i = l; i <= r; ++i)
                    A[vertices[i]] = false;

                if (l > 1)
                    A &= non_neighbors(vertices[vertices.size() - 1]);

                for (size_t i = l + 1; i < vertices.size() - 1; ++i) // i < vertices.size() - 2
                    A &= non_neighbors(vertices[i]);
                */

                auto A = neighbors(vertices[l]); // adjacent to p_l, non adjacent to p_{l+1}, ..., p_{k-1}, if l != 2: non adjacent to p_k
                for (size_t i = l + 1; i < k - 1; ++i)
                    A &= non_neighbors(vertices[i]);
                if (l != 1)
                    A &= non_neighbors(vertices[k - 1]);

                for (Vertex a : Graph::iterate(A)) {
                    vertices[l - 1] = a;
                    if (fill_outer(vertices, l - 1, vertices.size() - 1, graph, callback, neighbors, non_neighbors, valid_edge, valid_non_edge)) return true;
                }
            } else {
                /*
                auto C = graph.full_adjacency_row();
                for (size_t i = l; i <= r; ++i)
                    C[vertices[i]] = false;

                for (size_t i = l + 1; i < r - 1; ++i)
                    C &= non_neighbors(vertices[i]);

                auto A = C & neighbors(vertices[l]) & non_neighbors(vertices[r]);
                auto B = C & neighbors(vertices[r]) & non_neighbors(vertices[l]);
                */

                auto A = neighbors(vertices[l]); // adjacent to p_l, non adjacent to p_{l+1}, ..., p_r
                for (size_t i = l + 1; i <= r; ++i)
                    A &= non_neighbors(vertices[i]);

                auto B = neighbors(vertices[r]); // adjacent to p_r, non adjacent to p_l, ..., p_{r-1}
                for (size_t i = l; i < r; ++i)
                    B &= non_neighbors(vertices[i]);

                for (Vertex a : Graph::iterate(A)) {
                    for (Vertex b : Graph::iterate(B)) {
                        vertices[l - 1] = a;
                        vertices[r + 1] = b;
                        if (fill_outer(vertices, l - 1, r + 1, graph, callback, neighbors, non_neighbors, valid_edge, valid_non_edge)) return true;
                    }
                }
            }

            return false;
        }
    };

    template <>
    class FindNearImpl<3, false> {
    public:
        /**
         * Find paths of length 3 containing the vertices u and v. Each path is listed exactly once.
         *
         * @param graph
         * @param uv
         * @return
         */
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find_near(const Graph& graph, VertexPair uv, SubgraphCallback callback, F neighbors, G non_neighbors, H valid_edge, I valid_non_edge) {
            // assert(false);
            Graph::AdjRow Z(graph.size());
            auto [u, v] = uv;

            if (valid_edge(uv)) {
                Z = neighbors(v);
                Z &= non_neighbors(u);
                for (Vertex z : Graph::iterate(Z)) {
                    assert(valid_edge({u, v})); assert(valid_non_edge({u, z})); assert(valid_edge({v, z}));
                    if (callback(Subgraph{u, v, z})) return true;
                }

                Z = neighbors(u);
                Z &= non_neighbors(v);
                for (Vertex z : Graph::iterate(Z)) {
                    assert(valid_edge({z, u})); assert(valid_non_edge({z, v})); assert(valid_edge({u, v}));
                    if (callback(Subgraph{z, u, v})) return true;
                }
            } else if (valid_non_edge(uv)) {
                Z = neighbors(uv.u);
                Z &= neighbors(uv.v);
                for (Vertex z : Graph::iterate(Z)) {
                    assert(valid_edge({u, z})); assert(valid_non_edge({u, v})); assert(valid_edge({z, v}));
                    if (callback(Subgraph{u, z, v})) return true;
                }
            }
            return false;
        }
    };

    template <>
    class FindNearImpl<2, false> {
    public:
        /**
         * Find paths of length 2 containing vertices u and v. Each path is listed exactly once.
         *
         * @param uv
         * @return
         */
        template <typename SubgraphCallback, typename F, typename G, typename H, typename I>
        static bool find_near(const Graph& /*graph*/, VertexPair uv, SubgraphCallback callback, F /*neighbors*/, G /*non_neighbors*/, H valid_edge, I /*valid_non_edge*/) {
            assert(false);
            auto [u, v] = uv;

            if (valid_edge(uv)) {
                if (callback(Subgraph{u, v})) return true;
            }
            return false;
        }
    };
}

#endif //WEIGHTED_F_FREE_EDGE_EDITING_CENTERFINDNEARIMPL_H
