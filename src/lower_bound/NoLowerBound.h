//
// Created by jonas on 15.07.19.
//

#ifndef CONCEPT_NOLOWERBOUND_H
#define CONCEPT_NOLOWERBOUND_H

#include "../interfaces/LowerBoundI.h"

namespace LowerBound {
    class NoLowerBound : public LowerBoundI {
    private:
        class State : public StateI {
        public:
            std::unique_ptr<StateI> copy() override {
                return std::make_unique<State>(*this);
            }
        };

    public:
        explicit NoLowerBound(std::shared_ptr<FinderI> finder_ptr) : LowerBoundI(std::move(finder_ptr)) {};

        Cost result(StateI &/*state*/, Cost /*k*/) override {
            return 0;
        }

        std::unique_ptr<StateI> initialize(Cost /*k*/) override { return std::make_unique<State>(); }

        void before_mark_and_edit(StateI &/*state*/, VertexPair /*uv*/) override {}

        void after_mark_and_edit(StateI &/*state*/, VertexPair /*uv*/) override {}

        void before_mark(StateI &/*state*/, VertexPair /*uv*/) override {}

        void after_mark(StateI &/*state*/, VertexPair /*uv*/) override {}

        void before_edit(StateI &/*state*/, VertexPair /*uv*/) override {}

        void after_edit(StateI &/*state*/, VertexPair /*uv*/) override {}

        void after_unmark(StateI &/*state*/, VertexPair /*uv*/) override {}
    };
}


#endif //CONCEPT_NOLOWERBOUND_H
