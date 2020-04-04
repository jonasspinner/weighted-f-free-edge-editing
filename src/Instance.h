//
// Created by jonas on 17.07.19.
//

#ifndef CONCEPT_INSTANCE_H
#define CONCEPT_INSTANCE_H


#include <random>
#include "graph/VertexPairMap.h"

class Instance {
public:
    Graph graph;
    VertexPairMap<Cost> costs;

    std::string name;
    double multiplier;
    int permutation;

    Instance(std::string name_, Graph &&graph_, VertexPairMap<Cost> costs_, double multiplier_ = 1.0,
             int permutation_ = 0) :
            graph(std::move(graph_)), costs(std::move(costs_)), name(std::move(name_)), multiplier(multiplier_),
            permutation(permutation_) {}

    Instance(Instance &&other) : graph(std::move(other.graph)), costs(std::move(other.costs)),
                                 name(std::move(other.name)), multiplier(other.multiplier),
                                 permutation(other.permutation) {}

    Instance& operator=(Instance &&other) {
        using std::swap;
        swap(graph, other.graph);
        swap(costs, other.costs);
        swap(name, other.name);
        swap(multiplier, other.multiplier);
        swap(permutation, other.permutation);
        return *this;
    }

    Instance copy() const {
        return Instance(name, graph.copy(), costs, multiplier, permutation);
    }

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const Instance &instance) {
        using namespace YAML;
        out << BeginMap;
        out << Key << "name" << Value << instance.name;
        out << Key << "multiplier" << Value << instance.multiplier;
        out << Key << "permutation" << Value << instance.permutation;
        //out << Key << "graph" << Value << instance.graph;
        //out << Key << "costs" << Value << instance.costs;
        return out << EndMap;
    }
};


#endif //CONCEPT_INSTANCE_H
