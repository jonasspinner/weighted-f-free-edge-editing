#include "../src/lower_bound/LPRelaxation.h"
#include "../src/lower_bound/SortedGreedy.h"

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
                  const VertexPairMap<double> &values, Cost packing_lb, const VertexPairMap<bool> &covered_by_packing,
                  const std::vector<SubgraphT<Options::FSG::C4P4>> &packing,
                  const std::vector<VertexPair> &min_cost_vertex_pairs) {
    using namespace YAML;

    Emitter out;
    out << BeginDoc << BeginMap;
    out << Key << "relaxation_lb" << Value << relaxation_lb;
    out << Key << "packing_lb" << Value << packing_lb;
    output_matrix(out, "adjacency", adjacency);
    output_matrix(out, "costs", costs);
    output_matrix(out, "relaxation", values);
    output_matrix(out, "packing_cover", covered_by_packing);
    out << Key << "packing" << Value << packing;
    out << Key << "packing_min_cost_vertex_pairs" << Value << min_cost_vertex_pairs;
    out << EndMap << EndDoc;


    std::cout << out.c_str();
}


int main(int argc, char *argv[]) {
    namespace po = boost::program_options;
    constexpr Cost max_cost = std::numeric_limits<Cost>::max();


    Configuration config(Options::FSG::C4P4, 100, Options::SolverType::FPT, Options::Selector::FirstFound,
                         Options::LB::Trivial);
    config.input_path = "../data/bio/bio-nr-3-size-16.graph";

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("F", po::value<Options::FSG>(&config.forbidden_subgraphs)->default_value(config.forbidden_subgraphs))
            ("input", po::value<std::string>(&config.input_path)->default_value(config.input_path),
             "path to input instance")
            ("permutation", po::value<int>(&config.permutation)->default_value(config.permutation),
             "permutation of input instance")
            ("multiplier", po::value<double>(&config.multiplier)->default_value(config.multiplier),
             "multiplier for discretization of input instance");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }



    Instance instance = GraphIO::read_instance(config);
    VertexPairMap<bool> marked(instance.graph.size());

    auto edit_state = std::make_unique<EditState>(instance.graph.copy(), instance.costs);

    lower_bound::LPRelaxation<Options::FSG::C4P4> lp_algorithm(edit_state.get(), config.verbosity, config.timelimit);
    lp_algorithm.initialize(max_cost);
    auto relaxation_lb = lp_algorithm.calculate_lower_bound(max_cost);

    lower_bound::SortedGreedy<Options::FSG::C4P4> packing_algorithm(edit_state.get());
    packing_algorithm.initialize(max_cost);
    auto [packing_lb, covered_by_packing, packing, min_cost_vertex_pairs] = packing_algorithm.calculate_lower_bound_and_packing();

    VertexPairMap<bool> adjacency(instance.graph.size());
    for (VertexPair uv : instance.graph.edges()) {
        adjacency[uv] = true;
    }

    VertexPairMap<double> values(instance.graph.size());
    for (VertexPair uv : instance.graph.vertex_pairs()) {
        values[uv] = lp_algorithm.variable_edited_value(uv);
    }

    write_output(adjacency, instance.costs, relaxation_lb, values, packing_lb, covered_by_packing, packing, min_cost_vertex_pairs);


    return 0;
}

