//
// Created by jonas on 31.07.19.
//


#include "Center.h"

#include "FinderI.h"
#include "detail/CenterFindNearImpl.h"
#include "detail/CenterFindImpl.h"


namespace Finder {
    /**
     * Calls callback for all paths and cycles of the given length.
     *
     * @param graph
     * @param callback
     * @return
     */
    template <int length, bool with_cycles>
    bool Center<length, with_cycles>::find(const Graph& graph, SubgraphCallback callback) {
        auto cb = [&callback](Subgraph&& subgraph, Vertex) { return callback(std::move(subgraph)); };
        return detail::CenterFindImpl<length, with_cycles>::find(graph, cb,
                neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
    }

    /**
     * Calls callback for all paths and cycles of the given length. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
     *
     * @param graph
     * @param forbidden
     * @param callback
     * @return
     */
    template <int length, bool with_cycles>
    bool Center<length, with_cycles>::find(const Graph& graph, const Graph &forbidden, SubgraphCallback callback) {
        auto cb = [&callback](Subgraph&& subgraph, Vertex) { return callback(std::move(subgraph)); };
        return detail::CenterFindImpl<length, with_cycles>::find(graph, cb,
                neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    /**
     * Calls callback for all paths and cycles of the given length having both u and v as vertices.
     *
     * @param uv
     * @param graph
     * @param callback
     * @return
     */
    template <int length, bool with_cycles>
    bool Center<length, with_cycles>::find_near(VertexPair uv, const Graph& graph, SubgraphCallback callback) {
        // Note: Consider throwing an exception until efficient implementation is available.
        return detail::FindNearImpl<length, with_cycles>::find_near(graph, uv, callback,
                neighbors(graph), non_neighbors(graph), valid_edge(graph), valid_non_edge(graph));
    }

    /**
     * Calls callback for all paths and cycles of the given length having both u and v as vertices. Subgraphs sharing a vertex pair with the graph forbidden are ignored.
     *
     * @param uv
     * @param graph
     * @param forbidden
     * @param callback
     * @return
     */
    template <int length, bool with_cycles>
    bool Center<length, with_cycles>::find_near(VertexPair uv, const Graph& graph, const Graph &forbidden, SubgraphCallback callback) {
        // Note: Consider throwing an exception until efficient implementation is available.
        return detail::FindNearImpl<length, with_cycles>::find_near(graph, uv, callback,
                neighbors(graph, forbidden), non_neighbors(graph, forbidden), valid_edge(graph, forbidden), valid_non_edge(graph, forbidden));
    }

    template <int length, bool with_cycles>
    bool Center<length, with_cycles>::find_with_duplicates(const Graph &graph, const FinderI::SubgraphCallback &callback) {
        if constexpr (with_cycles) {
            return FinderI::find_with_duplicates(graph, callback);
        } else {
            return find(graph, callback);
        }
    }

    template <int length, bool with_cycles>
    bool Center<length, with_cycles>::find_with_duplicates(const Graph &graph, const Graph &forbidden,
                                                          const FinderI::SubgraphCallback &callback) {
        if constexpr (with_cycles) {
            return FinderI::find_with_duplicates(graph, forbidden, callback);
        } else {
            return find(graph, forbidden, callback);
        }
    }

    template <int length, bool with_cycles>
    bool Center<length, with_cycles>::for_all_conversionless_edits(const Subgraph &subgraph,
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
    Options::FSG Center<length, with_cycles>::forbidden_subgraphs() const {
        throw std::runtime_error("No forbidden subgraph available");
    }

    template <>
    Options::FSG Center<6, true>::forbidden_subgraphs() const {
        return Options::FSG::C6P6;
    }

    template <>
    Options::FSG Center<6, false>::forbidden_subgraphs() const {
        return Options::FSG::P6;
    }

    template <>
    Options::FSG Center<5, true>::forbidden_subgraphs() const {
        return Options::FSG::C5P5;
    }

    template <>
    Options::FSG Center<5, false>::forbidden_subgraphs() const {
        return Options::FSG::P5;
    }

    template <>
    Options::FSG Center<4, true>::forbidden_subgraphs() const {
        return Options::FSG::C4P4;
    }

    template <>
    Options::FSG Center<4, false>::forbidden_subgraphs() const {
        return Options::FSG::P4;
    }

    template <>
    Options::FSG Center<3, false>::forbidden_subgraphs() const {
        return Options::FSG::P3;
    }

    template <int length, bool with_cycles>
    std::string Center<length, with_cycles>::name() const {
        std::stringstream ss;
        ss << "CenterRec";
        if (with_cycles)
            ss << "C" << length;
        ss << "P" << length;
        return ss.str();
    }

    template <int length, bool with_cycles>
    void Center<length, with_cycles>::to_yaml(YAML::Emitter &out) const {
        using namespace YAML;
        out << BeginMap;
        out << Key << "name" << Value << "Center<" << length << ">";
        out << Key << "forbidden_subgraphs" << Value << "P" << length;
        if (with_cycles)
            out << "_C" << length;
        out << EndMap;
    }

    template class Center<6, true>;
    template class Center<5, true>;
    template class Center<4, true>;
    template class Center<6, false>;
    template class Center<5, false>;
    template class Center<4, false>;
    template class Center<3, false>;
}
