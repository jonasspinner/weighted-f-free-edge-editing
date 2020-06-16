#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LSSWZ_MWIS_SOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LSSWZ_MWIS_SOLVER_H


#include "LowerBoundI.h"
#include "../Instance.h"


namespace lower_bound {

    class LSSWZ_MWIS_Solver : public LowerBoundI {
    private:
        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;

        Configuration m_config;

        double time_limit = 10;
        bool disable_reduction = false;

        std::shared_ptr<FinderI> finder;
    public:
        LSSWZ_MWIS_Solver(const Instance &instance, const VertexPairMap<bool> &marked, Configuration config,
                          std::shared_ptr<FinderI> finder_ref) :
                m_graph(instance.graph), m_costs(instance.costs), m_marked(marked),
                m_config(std::move(config)), finder(std::move(finder_ref)) {}

        Cost calculate_lower_bound(Cost k) override;

    private:
        static std::optional<std::pair<Graph, std::vector<Cost>>>
        build_instance(FinderI &finder, const Graph &graph, const VertexPairMap<bool> &marked,
                       const VertexPairMap<Cost> &costs);
    };

}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_LSSWZ_MWIS_SOLVER_H
