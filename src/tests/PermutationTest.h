//
// Created by jonas on 20.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_PERMUTATIONTEST_H
#define WEIGHTED_F_FREE_EDGE_EDITING_PERMUTATIONTEST_H


class PermutationTest {
public:
    int seed;
    std::mt19937 gen;
    explicit PermutationTest(int seed_) : seed(seed_), gen(seed) {}

    void default_is_identity() {
        Permutation P(100);
        for (Vertex u = 0;  u < 100; ++u) {
            assert(u == P[u]);
        }
    }

    void reverse_is_correct() {
        Permutation P(100, seed);
        auto P_r = P.reverse();

        for (Vertex u = 0; u < 100; ++u) {
            assert(P_r[P[u]] == u);
        }
    }

    void graph_is_correct() {
        auto graph = random_graph(40, 100, gen);

        Permutation P(graph.size(), seed);
        auto P_r = P.reverse();

        std::vector<std::vector<Vertex>> adj(graph.size()), adj_p(graph.size());

        for (Vertex u : graph.vertices())
            for (Vertex v : graph.neighbors(u))
                adj[u].push_back(v);

        auto graph_p = P[graph];
        for (Vertex u : graph_p.vertices()) {
            for (Vertex v : graph_p.neighbors(u))
                adj_p[P_r[u]].push_back(P_r[v]);
            std::sort(adj_p[P_r[u]].begin(), adj_p[P_r[u]].end());
            }

        expect("Permutation of Graph is correct", adj, adj_p);
    }

    void run() {
        std::cout << "\nPermutationTest"
                     "\n---------------" << std::endl;
        default_is_identity();
        reverse_is_correct();
        graph_is_correct();
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_PERMUTATIONTEST_H
