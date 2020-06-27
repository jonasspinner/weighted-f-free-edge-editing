#ifndef WEIGHTED_F_FREE_EDGE_EDITING_LSSWZ_MWIS_SOLVER_H
#define WEIGHTED_F_FREE_EDGE_EDITING_LSSWZ_MWIS_SOLVER_H


#include "LowerBoundI.h"
#include "../Instance.h"
#include "../options.h"
#include "../Configuration.h"
#include "../forbidden_subgraphs/SubgraphC4P4.h"


namespace lower_bound {

    template<Options::FSG SetOfForbiddenSubgraphs>
    class LSSWZ_MWIS_Solver : public LowerBoundI {
    private:
        using Subgraph = SubgraphT<SetOfForbiddenSubgraphs>;
        using Finder = typename Subgraph::Finder;

        const Graph &m_graph;
        const VertexPairMap<Cost> &m_costs;
        const VertexPairMap<bool> &m_marked;

        Configuration m_config;

        double time_limit = 10;
        bool disable_reduction = false;

        Finder finder;
    public:
        LSSWZ_MWIS_Solver(const Instance &instance, const VertexPairMap<bool> &marked, Configuration config) :
                m_graph(instance.graph), m_costs(instance.costs), m_marked(marked),
                m_config(std::move(config)) {}

        Cost calculate_lower_bound(Cost k) override;

    private:
        static std::optional<std::pair<Graph, std::vector<Cost>>>
        build_instance(Finder &finder, const Graph &graph, const VertexPairMap<bool> &marked,
                       const VertexPairMap<Cost> &costs);
    };

}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_LSSWZ_MWIS_SOLVER_H
