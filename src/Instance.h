//
// Created by jonas on 17.07.19.
//

#ifndef CONCEPT_INSTANCE_H
#define CONCEPT_INSTANCE_H


#include "graph/VertexPairMap.h"

class Instance {
public:
    Graph graph;
    VertexPairMap<Cost> costs;

    Instance(Graph graph_, VertexPairMap<Cost> costs_) : graph(std::move(graph_)), costs(std::move(costs_)) {}
    Instance(Graph&& graph_, VertexPairMap<Cost>&& costs_) : graph(graph_), costs(costs_) {}
};


#endif //CONCEPT_INSTANCE_H
