//
// Created by jonas on 02.09.19.
//


#include "Endpoint.h"
#include "detail/EndpointFindImpl.h"


namespace Finder {
    template <int length, bool with_cycles>
    bool Endpoint<length, with_cycles>::find(const Graph& graph, SubgraphCallback callback) {
        auto cb = [&callback](Subgraph&& subgraph, Vertex) { return callback(std::move(subgraph)); };
        return detail::EndpointFindImpl<length, with_cycles>::find(graph, cb,
                                                                   neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
    }

    template <int length, bool with_cycles>
    bool Endpoint<length, with_cycles>::find(const Graph& graph, const Graph& forbidden, SubgraphCallback callback) {
        auto cb = [&callback](Subgraph&& subgraph, Vertex) { return callback(std::move(subgraph)); };
        return detail::EndpointFindImpl<length, with_cycles>::find(graph, cb,
                                                                   neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    template <int length, bool with_cycles>
    bool Endpoint<length, with_cycles>::find_near(VertexPair uv, const Graph& graph, SubgraphCallback callback) {
        // Note: Consider throwing an exception until efficient implementation is available.
        auto cb = [&](Subgraph&& subgraph, Vertex) {
            if (subgraph.contains(uv))
                return callback(std::move(subgraph));
            return false;
        };
        return detail::EndpointFindImpl<length, with_cycles>::find(graph, cb,
                                                                   neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
    }

    template <int length, bool with_cycles>
    bool Endpoint<length, with_cycles>::find_near(VertexPair uv, const Graph& graph, const Graph& forbidden, SubgraphCallback callback)  {
        // Note: Consider throwing an exception until efficient implementation is available.
        auto cb = [&](Subgraph&& subgraph, Vertex) {
            if (subgraph.contains(uv))
                return callback(std::move(subgraph));
            return false;
        };
        return detail::EndpointFindImpl<length, with_cycles>::find(graph, cb,
                                                                   neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }


    template <int length, bool with_cycles>
    bool Endpoint<length, with_cycles>::find_with_duplicates(const Graph &graph, const FinderI::SubgraphCallback &callback) {
        if constexpr (with_cycles) {
            return FinderI::find_with_duplicates(graph, callback);
        } else {
            return find(graph, callback);
        }
    }

    template <int length, bool with_cycles>
    bool Endpoint<length, with_cycles>::find_with_duplicates(const Graph &graph, const Graph &forbidden,
                                                          const FinderI::SubgraphCallback &callback) {
        if constexpr (with_cycles) {
            return FinderI::find_with_duplicates(graph, forbidden, callback);
        } else {
            return find(graph, forbidden, callback);
        }
    }

    template <int length, bool with_cycles>
    bool Endpoint<length, with_cycles>::for_all_conversionless_edits(const Subgraph &subgraph,
                                                                  const FinderI::VertexPairCallback &callback) const {
        if constexpr (with_cycles) {
            return FinderI::for_all_conversionless_edits(subgraph, callback);
        } else {
            for (auto uv : subgraph.vertexPairs()) {
                if (callback(uv))
                    return true;
            }
            return false;
        }
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
