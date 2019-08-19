//
// Created by jonas on 02.07.19.
//

#ifndef CONCEPT_VERTEXMAP_H
#define CONCEPT_VERTEXMAP_H


#include <iomanip>

#include "Graph.h"

template<typename T>
class VertexMap {
    using reference = typename std::vector<T>::reference;
    using const_reference = typename std::vector<T>::const_reference;

    Vertex n;
    std::vector<T> values;
public:
    explicit VertexMap(Vertex size) : n(size), values(size) {}

    const_reference operator[](Vertex vertex) const {
        assert(vertex < n);
        return values[vertex];
    }

    reference operator[](Vertex vertex) {
        assert(vertex < n);
        return values[vertex];
    }

    friend std::ostream &operator<<(std::ostream &os, const VertexMap &map) {
        os << "[";
        for (Vertex u = 0; u < map.n; ++u) {
            os << std::fixed << std::setw(8) << std::setprecision(4) << std::setfill(' ') << map[u] << " ";
        }
        return os << "]\n";
    }

    friend YAML::Emitter &operator<<(YAML::Emitter &out, const VertexMap &map) {
        out << YAML::BeginMap;
        out << YAML::Key << "size";
        out << YAML::Value << map.n;
        out << YAML::Key << "values";
        out << YAML::Value << map.values;
        out << YAML::EndMap;
        return out;
    }
};


#endif //CONCEPT_VERTEXMAP_H
