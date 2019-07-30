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

    Instance(Graph graph, VertexPairMap<Cost> costs) : graph(std::move(graph)), costs(std::move(costs)) {}
    Instance(Graph&& graph, VertexPairMap<Cost>&& costs) : graph(graph), costs(costs) {}
};


#endif //CONCEPT_INSTANCE_H
