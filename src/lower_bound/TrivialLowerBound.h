//
// Created by jonas on 15.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_TRIVIALLOWERBOUND_H
#define WEIGHTED_F_FREE_EDGE_EDITING_TRIVIALLOWERBOUND_H


#include "../interfaces/LowerBoundI.h"


namespace LowerBound {
    class TrivialLowerBound : public LowerBoundI {
    public:
        explicit TrivialLowerBound(std::shared_ptr<FinderI> finder_ptr) : LowerBoundI(std::move(finder_ptr)) {};

        Cost result(Cost /*k*/) override {
            return 0;
        }

        void push(Cost /*k*/) override {}

        void pop() override {}

        void before_mark_and_edit(VertexPair) override {}

        void after_mark_and_edit(VertexPair) override {}

        void before_mark(VertexPair) override {}

        void after_mark(VertexPair) override {}

        void before_edit(VertexPair) override {}

        void after_edit(VertexPair) override {}

        void after_unmark(VertexPair) override {}
    };
}


#endif //WEIGHTED_F_FREE_EDGE_EDITING_TRIVIALLOWERBOUND_H
