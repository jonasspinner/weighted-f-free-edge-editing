//
// Created by jonas on 27.02.20.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_INCREMENTBYMINCOST_H
#define WEIGHTED_F_FREE_EDGE_EDITING_INCREMENTBYMINCOST_H


class IncrementByMinCost : public SearchStrategyI {
public:
    Cost search_step() override { return -1; }
    void bound(Cost /*k*/, Cost /*lower_bound_k*/) override {}
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_INCREMENTBYMINCOST_H
