//
// Created by jonas on 03.07.19.
//

#ifndef CONCEPT_EXPLICITMWISSOLVER_H
#define CONCEPT_EXPLICITMWISSOLVER_H

#include <iostream>

#include "../graph/Graph.h"
#include "../graph/VertexMap.h"
#include "../graph/VertexPairMap.h"
#include "../finder/CenterC4P4.h"

namespace LowerBound {
    template<class Finder>
    class ExplicitMWISSolver : public LowerBoundI {
    private:
        struct MWISInstance {
            Graph graph;
            VertexMap<float> weights;

            MWISInstance(Graph &&graph, VertexMap<float> &&weights) : graph(std::move(graph)),
                                                                      weights(std::move(weights)) {}
        };

        std::shared_ptr<FinderI> finder;
        const Graph &graph;
        const VertexPairMap<float> &weights;

    public:
        ExplicitMWISSolver(const Graph &graph, const VertexPairMap<float> &weights, std::shared_ptr<FinderI> finder)
                : finder(std::move(finder)), graph(graph),
                  weights(weights) {};

        size_t update(VertexPair pair) {
            std::vector<Subgraph> subgraphs;
            finder->find([&](const auto &subgraph) {
                subgraphs.push_back(subgraph);
                return false;
            });

            auto instance = construct(graph, weights, subgraphs);

            std::cout << instance.graph;
            std::cout << instance.weights;

            std::vector<Vertex> bound;
            VertexMap<bool> marked(instance.graph);
            for (Vertex i = 0; i < instance.graph.size(); ++i) {
                if (!marked[i]) {
                    bound.push_back(i);
                    instance.graph.for_neighbors_of(i, [&](auto j) {
                        marked[j] = true;
                        return false;
                    });
                }
            }
            std::cout << "greedy lower bound: " << bound.size() << "\n";
            for (auto i : bound) {
                std::cout << i << ":" << instance.weights[i] << " ";
            }
            std::cout << "\n";

            return 0;
        }

    private:
        static MWISInstance construct(const Graph &graph, const VertexPairMap<float> &vertex_pair_weights,
                                      const std::vector<Subgraph> &subgraphs) {
            VertexPairMap<std::vector<Vertex>> edges_to_sg(graph.size());
            for (size_t i = 0; i < subgraphs.size(); ++i) {
                for (VertexPair uv : subgraphs[i].vertexPairs()) {
                    edges_to_sg[uv].push_back(i);
                }
            }

            Graph MWISGraph(subgraphs.size());
            VertexMap<float> weights(subgraphs.size());

            for (size_t i = 0; i < subgraphs.size(); ++i) {
                weights[i] = min_vertex_pair_weight(vertex_pair_weights, subgraphs[i]);
            }

            graph.for_all_vertex_pairs([&](const auto &edge) {
                const auto &sgs = edges_to_sg[edge];
                for (int i = 0; i < sgs.size(); ++i) {
                    for (int j = i + 1; j < sgs.size(); ++j) {
                        auto u = sgs[i], v = sgs[j];
                        MWISGraph.set_edge({u, v});
                    }
                }
                return false;
            });

            return MWISInstance(std::move(MWISGraph), std::move(weights));
        }

        static float
        min_vertex_pair_weight(const VertexPairMap<float> &vertex_pair_weights, const Subgraph &subgraph) {
            float min_value = std::numeric_limits<float>::max();
            for (VertexPair uv : subgraph.vertexPairs()) {
                auto value = vertex_pair_weights[uv];
                if (value < min_value) min_value = value;
            }
            return min_value;
        }
    };
}


#endif //CONCEPT_EXPLICITMWISSOLVER_H
