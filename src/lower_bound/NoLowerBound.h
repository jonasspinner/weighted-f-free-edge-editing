//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_NOLOWERBOUND_H
#define CONCEPT_NOLOWERBOUND_H

#include "../interfaces/LowerBoundI.h"

namespace LowerBound {
    class NoLowerBound : public LowerBoundI {
    private:

    public:
        explicit NoLowerBound(std::shared_ptr<FinderI> finder_ptr) : LowerBoundI(std::move(finder_ptr)) {};

        Cost result(Cost /*k*/) override {
            return 0;
        }

        void push(Cost k) override {}

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


#endif //CONCEPT_NOLOWERBOUND_H
