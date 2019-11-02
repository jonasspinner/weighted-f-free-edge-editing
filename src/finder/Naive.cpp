//
// Created by jonas on 02.11.19.
//

#include "Naive.h"
#include "detail/NaiveFindImpl.h"


namespace Finder {
    template <int length, bool with_cycles>
    bool Naive<length, with_cycles>::find(SubgraphCallback callback) {
        return detail::NaiveFindImpl<length, with_cycles>::find(graph, callback, valid_edge(graph), valid_non_edge(graph));
    }

    template <int length, bool with_cycles>
    bool Naive<length, with_cycles>::find(const Graph& forbidden, SubgraphCallback callback) {
        return detail::NaiveFindImpl<length, with_cycles>::find(graph, callback, valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    template <int length, bool with_cycles>
    bool Naive<length, with_cycles>::find_near(VertexPair uv, SubgraphCallback callback) {
        auto cb = [&](Subgraph&& subgraph) {
            if (subgraph.contains(uv))
                return callback(std::move(subgraph));
            return false;
        };
        return detail::NaiveFindImpl<length, with_cycles>::find(graph, cb, valid_edge(graph), valid_non_edge(graph));
    }

    template <int length, bool with_cycles>
    bool Naive<length, with_cycles>::find_near(VertexPair uv, const Graph& forbidden, SubgraphCallback callback)  {
        auto cb = [&](Subgraph&& subgraph) {
            if (subgraph.contains(uv))
                return callback(std::move(subgraph));
            return false;
        };
        return detail::NaiveFindImpl<length, with_cycles>::find(graph, cb, valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    template <int length, bool with_cycles>
    Options::FSG Naive<length, with_cycles>::forbidden_subgraphs() const {
        throw std::runtime_error("No forbidden subgraph available");
    }

    template <>
    Options::FSG Naive<6, true>::forbidden_subgraphs() const {
        return Options::FSG::C6P6;
    }

    template <>
    Options::FSG Naive<6, false>::forbidden_subgraphs() const {
        return Options::FSG::P6;
    }

    template <>
    Options::FSG Naive<5, true>::forbidden_subgraphs() const {
        return Options::FSG::C5P5;
    }

    template <>
    Options::FSG Naive<5, false>::forbidden_subgraphs() const {
        return Options::FSG::P5;
    }

    template <>
    Options::FSG Naive<4, true>::forbidden_subgraphs() const {
        return Options::FSG::C4P4;
    }

    template <>
    Options::FSG Naive<4, false>::forbidden_subgraphs() const {
        return Options::FSG::P4;
    }

    template <>
    Options::FSG Naive<3, false>::forbidden_subgraphs() const {
        return Options::FSG::P3;
    }

    template <int length, bool with_cycles>
    std::string Naive<length, with_cycles>::name() const {
        std::stringstream ss;
        ss << "NaiveRec";
        if (with_cycles)
            ss << "C" << length;
        ss << "P" << length;
        return ss.str();
    }

    template <int length, bool with_cycles>
    void Naive<length, with_cycles>::to_yaml(YAML::Emitter &out) const {
        using namespace YAML;
        out << BeginMap;
        out << Key << "name" << Value << "Naive<" << length << ">";
        out << Key << "forbidden_subgraphs" << Value << "P" << length;
        if (with_cycles)
            out << "_C" << length;
        out << EndMap;
    }

    template class Naive<6, true>;
    template class Naive<5, true>;
    template class Naive<4, true>;
    template class Naive<6, false>;
    template class Naive<5, false>;
    template class Naive<4, false>;
    template class Naive<3, false>;
}
