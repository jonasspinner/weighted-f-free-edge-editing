//
// Created by jonas on 09.03.20.
//


#include "../src/lower_bound/LPRelaxation.h"
#include "../src/lower_bound/SortedGreedy.h"

#include "../src/finder/finder_utils.h"
#include "../src/graph/GraphIO.h"


VertexPairMap<double> to_double(const VertexPairMap<bool> &map) {
    VertexPairMap<double> new_map(map.size());
    for (auto uv : Graph::VertexPairs(map.size())) {
        new_map[uv] = static_cast<double>(map[uv]);
    }
    return new_map;
}

template<typename T>
void output_matrix(YAML::Emitter &out, const std::string &key, const VertexPairMap<T> &map) {
    using namespace YAML;
    out << Key << key << Value << BeginMap;
    out << Key << "size" << Value << map.size();
    out << Key << "values" << Value << BeginSeq;
    for (Vertex u : Graph::Vertices(map.size())) {
        out << Flow << BeginSeq;
        for (Vertex v : Graph::Vertices(map.size())) {
            if (u == v) {
                out << T();
            } else {
                out << map[{u, v}];
            }
        }
        out << EndSeq;
    }
    out << EndSeq << EndMap;
}

void output_matrix(YAML::Emitter &out, const std::string &key, const VertexPairMap<bool> &map) {
    auto new_map = to_double(map);
    output_matrix(out, key, new_map);
}


void write_output(const VertexPairMap<bool> &adjacency, const VertexPairMap<Cost> &costs, Cost relaxation_lb,
                  const VertexPairMap<double> &values, Cost packing_lb, const VertexPairMap<bool> &used_by_packing) {
    using namespace YAML;

    Emitter out;
    out << BeginDoc << BeginMap;
    out << Key << "relaxation_lb" << Value << relaxation_lb;
    out << Key << "packing_lb" << Value << packing_lb;
    output_matrix(out, "adjacency", adjacency);
    output_matrix(out, "costs", costs);
    output_matrix(out, "relaxation", values);
    output_matrix(out, "packing", used_by_packing);
    out << EndMap << EndDoc;


    std::cout << out.c_str();
}


int main(int argc, char *argv[]) {

    constexpr Cost max_cost = std::numeric_limits<Cost>::max();


    Configuration config(Options::FSG::C4P4, 100, Options::SolverType::FPT, Options::Selector::FirstFound,
                         Options::LB::Trivial);
    config.input_path = "../data/bio/bio-nr-3-size-16.graph";

    Instance instance = GraphIO::read_instance(config);
    VertexPairMap<bool> marked(instance.graph.size());


    std::shared_ptr<FinderI> finder = Finder::make(config.forbidden_subgraphs);

    lower_bound::LPRelaxation lp_algorithm(instance, marked, config, finder);
    lp_algorithm.initialize(max_cost);
    auto relaxation_lb = lp_algorithm.calculate_lower_bound(max_cost);

    lower_bound::SortedGreedy packing_algorithm(instance, marked, finder);
    packing_algorithm.initialize(max_cost);
    auto packing_lb = packing_algorithm.calculate_lower_bound(max_cost);

    //std::cout << relaxation_lb << "\n";
    //std::cout << packing_lb << "\n";

    VertexPairMap<bool> adjacency(instance.graph.size());
    for (VertexPair uv : instance.graph.edges()) {
        adjacency[uv] = true;
    }

    VertexPairMap<double> values(instance.graph.size());
    for (VertexPair uv : instance.graph.vertexPairs()) {
        values[uv] = lp_algorithm.variable_edited_value(uv);
    }

    write_output(adjacency, instance.costs, relaxation_lb, values, packing_lb, packing_algorithm.used_in_bound());


    return 0;
}

