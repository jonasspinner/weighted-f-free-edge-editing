//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_CONSUMERI_H
#define CONCEPT_CONSUMERI_H


#include "../graph/Graph.h"

#include "../interfaces/FinderI.h"


class ConsumerI {
protected:
    std::shared_ptr<FinderI> finder;

public:
    explicit ConsumerI(std::shared_ptr<FinderI> finder_ptr) : finder(std::move(finder_ptr)) {}

    virtual ~ConsumerI() = default;

    // virtual std::unique_ptr<StateI> initialize(Cost k) = 0;

protected:
    virtual void push(Cost k) = 0;

    virtual void pop() = 0;

public:
    // all
    virtual void before_mark_and_edit(VertexPair uv) = 0;

    // all
    virtual void after_mark_and_edit(VertexPair uv) = 0;

    // all
    virtual void before_mark(VertexPair uv) = 0;

    // subgraph_stats
    virtual void after_mark(VertexPair uv) = 0;

    // subgraph_stats
    virtual void before_edit(VertexPair uv) = 0;

    // subgraph_stats
    virtual void after_edit(VertexPair uv) = 0;

    // subgraph_stats
    virtual void after_unmark(VertexPair uv) = 0;

    friend class Level;
};

#endif //CONCEPT_CONSUMERI_H
