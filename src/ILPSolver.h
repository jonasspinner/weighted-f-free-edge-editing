//
// Created by jonas on 19.08.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ILPSOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ILPSOLVER_H


#include <gurobi_c++.h>

#include "finder/finder_utils.h"
#include "../src/Solver.h"
#include "../src/Instance.h"
#include "Configuration.h"


class ILPSolver : public Solver {
private:
    class FSGCallback : public GRBCallback {
    private:
        Graph m_graph;
        std::unique_ptr<FinderI> m_finder;
        const VertexPairMap<GRBVar> &m_vars;
    public:
        FSGCallback(Options::FSG fsg, Graph graph, const VertexPairMap<GRBVar> &vars) : m_graph(std::move(graph)),
                                                                                        m_vars(vars) {
            m_finder = Finder::make(fsg, m_graph);
        };

    protected:
        void callback() override {
            if (where == GRB_CB_MIPSOL) {

                VertexPairMap<bool> edited(m_graph.size());
                std::vector<VertexPair> edits;

                // List edits of current solution. Apply them to graph.
                for (VertexPair uv : m_graph.vertexPairs())
                    if (getSolution(m_vars[uv]) > 0.5) {
                        edits.push_back(uv);
                        m_graph.toggleEdge(uv);
                        edited[uv] = true;
                    }

                // Find forbidden subgraphs in current solution and add additional constraints.
                m_finder->find([&](Subgraph &&subgraph) {
                    GRBLinExpr pairs;
                    for (VertexPair uv : subgraph.vertexPairs()) {
                        if (edited[uv]) {
                            pairs += (1 - m_vars[uv]);
                        } else {
                            pairs += m_vars[uv];
                        }
                    }

                    addLazy(pairs >= 1);
                    return false;
                });

                // Undo edits.
                for (VertexPair uv : edits)
                    m_graph.toggleEdge(uv);

            }
        }
    };

    Options::FSG m_fsg;

public:
    explicit ILPSolver(Options::FSG fsg) : m_fsg(fsg) {}

    std::vector<Solution> solve(Instance instance) override {
        auto &[_, graph, costs] = instance;

        try {
            GRBEnv env = GRBEnv();
            GRBModel model = GRBModel(env);

            model.set(GRB_IntParam_LazyConstraints, 1);


            // Add model variables: x_uv == 1  <=>  uv is edited
            VertexPairMap<GRBVar> vars(graph.size());
            for (VertexPair uv : graph.vertexPairs())
                vars[uv] = model.addVar(0.0, 1.0, 0, GRB_BINARY);


            // Set objective. When uv is edited the objective is increased by costs[uv]
            GRBLinExpr obj;
            for (VertexPair uv : graph.vertexPairs())
                obj += costs[uv] * vars[uv];

            model.setObjective(obj, GRB_MINIMIZE);


            // Register callback
            FSGCallback cb(m_fsg, graph, vars);
            model.setCallback(&cb);


            // Optimize model
            model.optimize();


            // Read solution
            std::vector<VertexPair> edits;
            for (VertexPair uv : graph.vertexPairs())
                if (is_set(vars, uv))
                    edits.push_back(uv);

            return {Solution(instance, edits)};

        } catch (GRBException &e) {
            std::cerr << "Error " << e.getErrorCode() << " " << e.getMessage();
            abort();
        }
    }

private:
    static bool is_set(const VertexPairMap<GRBVar> &vars, VertexPair uv) {
        double value = vars[uv].get(GRB_DoubleAttr_X);
        assert(0.01 < value || value < 0.99);
        return value > 0.5;
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ILPSOLVER_H
