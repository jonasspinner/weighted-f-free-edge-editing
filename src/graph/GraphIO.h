//
// Created by jonas on 02.07.19.
//

#ifndef CONCEPT_GRAPHIO_H
#define CONCEPT_GRAPHIO_H

#include <string>
#include <fstream>

#include "Graph.h"
#include "VertexPairMap.h"
#include "../Instance.h"
#include "../Permutation.h"

class GraphIO {
public:
    struct Format {
        struct Metis {
        };
    };

    /**
     * Reads a graph from path.
     *
     * Assumes that the file is in metis format with fmt == 1 and that the edges are from the upper triangular adjacency
     * matrix of a fully connected graph.
     *
     * The instance is constructed as follows. Editing costs c(uv) = ceil(abs(s(uv)) * multiplier). If permutation is
     * not 0, then the instance is permuted with the given permutation as seed.
     *
     * @param path The path to the instance file.
     * @param multiplier A factor for multiplying the given similarity scores.
     * @param permutation
     * @return
     */
    static Instance read_instance(const std::string &path, double multiplier = 1.0, int permutation = 0) {
        using Real = double;

        std::ifstream file(path);

        if (!file){
            std::stringstream ss;
            ss << "could not open file " << path;
            throw std::runtime_error(ss.str());
        }

        auto next_non_comment_line = [](std::ifstream &f, std::string &line) {
            do {
                std::getline(f, line);
            } while (!line.empty() && line[0] == '%');
        };

        Vertex n, m, fmt;
        std::string line;
        next_non_comment_line(file, line);
        std::stringstream(line) >> n >> m >> fmt;


        Graph G(n);
        Cost undefined_cost = std::numeric_limits<Cost>::max();
        VertexPairMap<Cost> edit_costs(G.size(), undefined_cost);

        if (fmt == 0) {
            edit_costs = VertexPairMap<Cost>(G.size(), 1);
            for (Vertex u = 1; u <= n; ++u) {
                next_non_comment_line(file, line);
                std::stringstream ss(line);
                Vertex v;
                while (ss >> v) {
                    if (u <= v) {
                        if (u == v)
                            throw std::runtime_error("self loops are not allowed");
                        G.setEdge({u - 1, v - 1});
                    }
                }
            }
        } else if (fmt == 1) {
            for (Vertex u = 1; u <= n; ++u) {
                next_non_comment_line(file, line);
                std::stringstream ss(line);
                Vertex v;
                Real weight;
                while (ss >> v >> weight) {
                    if (u <= v) {
                        if (u == v)
                            throw std::runtime_error("self loops are not allowed");

                        VertexPair edge(u - 1, v - 1);
                        Cost edge_cost = static_cast<Cost>(std::ceil(std::abs(weight) * multiplier));

                        if (edit_costs[edge] != undefined_cost && edit_costs[edge] != edge_cost)
                            throw std::runtime_error("weights are not symmetric");

                        edit_costs[edge] = edge_cost;
                        if (weight >= 0) G.setEdge(edge);
                        else G.clearEdge(edge);
                    }
                }
            }
        } else {
            throw std::runtime_error("fmt not supported");
        }

        for (VertexPair uv : G.vertexPairs())
            if (edit_costs[uv] == undefined_cost)
                throw std::runtime_error("undefined editing cost");

        Instance instance(path, G, edit_costs, multiplier);
        Permutation P(G.size(), permutation);
        instance = P[instance];

        return instance;
    }

    template<typename Weight>
    static void write_instance(const std::string &path, const Graph &graph, const VertexPairMap<Weight> &weights) {
        write_instance(path, graph, weights, Format::Metis());
    }

    template<typename Weight>
    static void
    write_instance(const std::string &path, const Graph &graph, const VertexPairMap<Weight> &weights, Format::Metis) {
        std::ofstream file(path);

        if (!file) throw std::runtime_error("could not open file");

        Vertex n = graph.size();

        file << n << " " << n * (n - 1) / 2 << " " << 1 << "\n";
        for (Vertex u : graph.vertices()) {
            for (Vertex v : graph.vertices()) {
                if (u >= v) continue;
                file << (v + 1) << " " << (graph.hasEdge({u, v}) ? 1 : -1) * weights[{u, v}]
                     << " ";
            }
            file << "\n";
        }
    }
};


#endif //CONCEPT_GRAPHIO_H
