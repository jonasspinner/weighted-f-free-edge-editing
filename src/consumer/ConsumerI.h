//
// Created by jonas on 15.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_CONSUMERI_H
#define WEIGHTED_F_FREE_EDGE_EDITING_CONSUMERI_H


#include "../graph/Graph.h"

#include "../finder/FinderI.h"


class ConsumerI {
public:
    virtual ~ConsumerI() = default;

    virtual void initialize(Cost /*k*/) {};

    virtual void push_state(Cost /*k*/) {}; // push next_state = copy(state)

    virtual void pop_state() {}; // pop state

    virtual void before_mark(VertexPair /*uv*/) {};             // all

    virtual void after_mark(VertexPair /*uv*/) {};              // subgraph_stats

    virtual void before_edit(VertexPair /*uv*/) {};             // subgraph_stats

    virtual void after_edit(VertexPair /*uv*/) {};              // subgraph_stats

    virtual void before_unedit(VertexPair /*uv*/) {};           // subgraph_stats

    virtual void after_unedit(VertexPair /*uv*/) {};            // subgraph_stats

    virtual void after_unmark(VertexPair /*uv*/) {};            // subgraph_stats

    friend class Editor;
};

#endif //WEIGHTED_F_FREE_EDGE_EDITING_CONSUMERI_H
