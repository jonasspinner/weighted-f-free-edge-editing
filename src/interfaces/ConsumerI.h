//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_CONSUMERI_H
#define CONCEPT_CONSUMERI_H


#include "../graph/Graph.h"

#include "../interfaces/FinderI.h"


class StateI {
public:
    virtual ~StateI() = default;
    virtual std::unique_ptr<StateI> copy() = 0;
};


class ConsumerI {
protected:
    std::shared_ptr<FinderI> finder;

public:
    explicit ConsumerI(std::shared_ptr<FinderI> finder_ptr) : finder(std::move(finder_ptr)) {}

    virtual ~ConsumerI() = default;

    virtual std::unique_ptr<StateI> initialize(Cost k) = 0;

    // all
    virtual void before_mark_and_edit(StateI &state, VertexPair uv) = 0;

    // all
    virtual void after_mark_and_edit(StateI &state, VertexPair uv) = 0;

    // all
    virtual void before_mark(StateI &state, VertexPair uv) = 0;

    // subgraph_stats
    virtual void after_mark(StateI &state, VertexPair uv) = 0;

    // subgraph_stats
    virtual void before_edit(StateI &state, VertexPair uv) = 0;

    // subgraph_stats
    virtual void after_edit(StateI &state, VertexPair uv) = 0;

    // subgraph_stats
    virtual void after_unmark(StateI &state, VertexPair uv) = 0;

};

#endif //CONCEPT_CONSUMERI_H
