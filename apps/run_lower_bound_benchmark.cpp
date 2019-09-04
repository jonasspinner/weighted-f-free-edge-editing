//
// Created by jonas on 04.09.19.
//


#include <iostream>
#include <chrono>

#include "../src/graph/GraphIO.h"
#include "../src/Configuration.h"

#include "../src/lower_bound/utils.h"
#include "../src/finder/utils.h"


int main() {
    const std::vector<std::string> paths {
            "./data/bio/bio-nr-3-size-16.metis",
            "./data/bio/bio-nr-4-size-39.metis",
            "./data/bio/bio-nr-11-size-22.metis",
            "./data/bio/bio-nr-277-size-222.metis",
            "./data/karate.graph",
            "./data/lesmis.graph",
            "./data/dolphins.graph",
            "./data/grass_web.metis.graph"};


    using Options::LB;
    std::vector<Options::LB> lower_bounds = {LB::Greedy, LB::No, LB::LocalSearch, LB::LinearProgram};


    auto instance = GraphIO::read_graph(paths[0], 100);
    VertexPairMap<bool> marked(instance.graph.size());

    std::shared_ptr<FinderI> finder = Finder::make(Options::FSG::P4C4, instance.graph);
    
    constexpr Cost max_k = std::numeric_limits<Cost>::max();

    for (auto lower_bound : lower_bounds) {
        auto lb = LowerBound::make(lower_bound, finder, instance, marked);
        

        auto t1 = std::chrono::steady_clock::now();
        lb->initialize();
        auto t2 = std::chrono::steady_clock::now();
        auto value = lb->result(max_k);
        auto t3 = std::chrono::steady_clock::now();
        
        using namespace YAML;
        Emitter out;
        out << BeginDoc << BeginMap;
        out << Key << "lower_bound_name" << Value << lower_bound;
        out << Key << "instance" << Value << instance.name;
        out << Key << "value" << Value << value;
        out << Key << "initialization_time" << Value << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() << Comment("ns");
        out << Key << "result_time" << Value << std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count() << Comment("ns");
        out << Key << "complete_time" << Value << std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t1).count() << Comment("ns");
        out << EndMap << EndDoc;
        std::cout << out.c_str();
    }

    return 0;
}

