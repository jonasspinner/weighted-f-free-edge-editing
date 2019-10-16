//
// Created by jonas on 02.09.19.
//


#include "Endpoint.h"
#include "detail/EndpointFindImpl.h"


namespace Finder {
    template <int length, bool with_cycles>
    bool Endpoint<length, with_cycles>::find(SubgraphCallback callback) {
        auto cb = [&callback](Subgraph&& subgraph, Vertex) { return callback(std::move(subgraph)); };
        return detail::EndpointFindImpl<length, with_cycles>::find(graph, cb,
                                                                   neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
    }

    template <int length, bool with_cycles>
    bool Endpoint<length, with_cycles>::find(const Graph& forbidden, SubgraphCallback callback) {
        auto cb = [&callback](Subgraph&& subgraph, Vertex) { return callback(std::move(subgraph)); };
        return detail::EndpointFindImpl<length, with_cycles>::find(graph, cb,
                                                                   neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    template <int length, bool with_cycles>
    bool Endpoint<length, with_cycles>::find_near(VertexPair uv, SubgraphCallback callback) {
        auto cb = [&](Subgraph&& subgraph, Vertex) {
            if (subgraph.contains(uv))
                return callback(std::move(subgraph));
            return false;
        };
        return detail::EndpointFindImpl<length, with_cycles>::find(graph, cb,
                                                                   neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
    }

    template <int length, bool with_cycles>
    bool Endpoint<length, with_cycles>::find_near(VertexPair uv, const Graph& forbidden, SubgraphCallback callback)  {
        auto cb = [&](Subgraph&& subgraph, Vertex) {
            if (subgraph.contains(uv))
                return callback(std::move(subgraph));
            return false;
        };
        return detail::EndpointFindImpl<length, with_cycles>::find(graph, cb,
                                                                   neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    template <int length, bool with_cycles>
    Options::FSG Endpoint<length, with_cycles>::forbidden_subgraphs() const {
        throw std::runtime_error("No forbidden subgraph available");
    }

    template <>
    Options::FSG Endpoint<6, true>::forbidden_subgraphs() const {
        return Options::FSG::C6P6;
    }

    template <>
    Options::FSG Endpoint<6, false>::forbidden_subgraphs() const {
        return Options::FSG::P6;
    }

    template <>
    Options::FSG Endpoint<5, true>::forbidden_subgraphs() const {
        return Options::FSG::C5P5;
    }

    template <>
    Options::FSG Endpoint<5, false>::forbidden_subgraphs() const {
        return Options::FSG::P5;
    }

    template <>
    Options::FSG Endpoint<4, true>::forbidden_subgraphs() const {
        return Options::FSG::C4P4;
    }

    template <>
    Options::FSG Endpoint<4, false>::forbidden_subgraphs() const {
        return Options::FSG::P4;
    }
    
    template <>
    Options::FSG Endpoint<3, false>::forbidden_subgraphs() const {
        return Options::FSG::P3;
    }
    
    template <int length, bool with_cycles>
    std::string Endpoint<length, with_cycles>::name() const {
        std::stringstream ss;
        ss << "EndpointRec";
        if (with_cycles)
            ss << "C" << length;
        ss << "P" << length;
        return ss.str();
    }

    template <int length, bool with_cycles>
    void Endpoint<length, with_cycles>::to_yaml(YAML::Emitter &out) const {
        using namespace YAML;
        out << BeginMap;
        out << Key << "name" << Value << "Endpoint<" << length << ">";
        out << Key << "forbidden_subgraphs" << Value << "P" << length;
        if (with_cycles)
            out << "_C" << length;
        out << EndMap;
    }

    template class Endpoint<6, true>;
    template class Endpoint<5, true>;
    template class Endpoint<4, true>;
    template class Endpoint<6, false>;
    template class Endpoint<5, false>;
    template class Endpoint<4, false>;
    template class Endpoint<3, false>;
}
