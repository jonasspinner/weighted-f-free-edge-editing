//
// Created by jonas on 17.08.19.
//

#include <sstream>
#include <utility>
#include <gurobi_c++.h>

#include "graph/GraphIO.h"
#include "finder/CenterC4P4.h"


class FSGCallback : public GRBCallback {
private:
    Graph graph;
    Finder::CenterC4P4 finder;
    const VertexPairMap<GRBVar> &vars;
public:
    FSGCallback(Graph graph_, const VertexPairMap<GRBVar> &vars_) : graph(std::move(graph_)), finder(graph), vars(vars_) {};

protected:
    void callback() override {
        if (where == GRB_CB_MIPSOL) {
            VertexPairMap<bool> edited(graph.size());
            std::vector<VertexPair> edits;
            for (VertexPair uv : graph.vertexPairs())
                if (getSolution(vars[uv]) > 0.5) {
                    edits.push_back(uv);
                    graph.toggle_edge(uv);
                    edited[uv] = true;
                }

            finder.find([&](Subgraph &&subgraph) {
                GRBLinExpr pairs;
                for (VertexPair uv : subgraph.vertexPairs()) {
                    if (edited[uv]) {
                        pairs += (1 - vars[uv]);
                    } else {
                        pairs += vars[uv];
                    }
                }

                addLazy(pairs >= 1);
                return false;
            });

            for (VertexPair uv : edits)
                graph.toggle_edge(uv);
        }
    }
};

int main() {
    const std::vector<std::string> paths {
            "../data/cost_matrix_component_nr_3_size_16_cutoff_10.0.metis",
            "../data/cost_matrix_component_nr_4_size_39_cutoff_10.0.metis",
            "../data/cost_matrix_component_nr_11_size_22_cutoff_10.0.metis",
            "../data/cost_matrix_component_nr_277_size_222_cutoff_10.0.metis",
            "./data/karate.graph",
            "./data/lesmis.graph",
            "./data/dolphins.graph",
            "./data/grass_web.metis.graph"};
    auto instance = GraphIO::read_graph(paths[6], 100);

    auto &graph = instance.graph;
    auto &costs = instance.costs;

    //for (VertexPair uv : graph.vertexPairs())
    //    costs[uv] = 1;
    // std::cout << costs;

    Finder::CenterC4P4 finder(graph);


    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    model.set(GRB_IntParam_LazyConstraints, 1);

    VertexPairMap<GRBVar> vars(graph.size());
    for (VertexPair uv : graph.vertexPairs()) {
        std::stringstream name;
        name << "edit " << uv;
        vars[uv] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name.str());
    }

    GRBLinExpr obj;
    for (VertexPair uv : graph.vertexPairs())
        obj += costs[uv] * vars[uv];

    model.setObjective(obj, GRB_MINIMIZE);

    FSGCallback cb(graph, vars);
    model.setCallback(&cb);

    // Optimize model
    model.optimize();

    for (VertexPair uv : graph.vertexPairs()) {
        if (vars[uv].get(GRB_DoubleAttr_X) > 0.5) {
            std::cout << vars[uv].get(GRB_StringAttr_VarName) << " ";
            graph.toggle_edge(uv);
        }
    }
    std::cout << "\n";

    finder.find([&](Subgraph &&subgraph) {
        std::cerr << subgraph << "\n";
        return false;
    });


    std::cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;


    return 0;
}