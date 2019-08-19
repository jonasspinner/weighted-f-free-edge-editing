//
// Created by jonas on 17.07.19.
//

#ifndef CONCEPT_INSTANCE_H
#define CONCEPT_INSTANCE_H


#include <random>
#include "graph/VertexPairMap.h"

class Instance {
public:
    std::string name;
    Graph graph;
    VertexPairMap<Cost> costs;

    Instance(std::string name_, Graph graph_, VertexPairMap<Cost> costs_) : name(std::move(name_)), graph(std::move(graph_)), costs(std::move(costs_)) {}


    friend YAML::Emitter &operator<<(YAML::Emitter &out, const Instance &instance) {
        out << YAML::BeginMap;
        out << YAML::Key << "name";
        out << YAML::Value << instance.name;
        out << YAML::Key << "graph";
        out << YAML::Value << instance.graph;
        out << YAML::Key << "costs";
        out << YAML::Value << instance.costs;
        return out << YAML::EndMap;
    }
};


#endif //CONCEPT_INSTANCE_H
