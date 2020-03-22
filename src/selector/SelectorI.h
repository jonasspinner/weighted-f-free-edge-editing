//
// Created by jonas on 17.07.19.
//

#ifndef WEIGHTED_F_FREE_EDGE_EDITING_SELECTORI_H
#define WEIGHTED_F_FREE_EDGE_EDITING_SELECTORI_H

#include "../consumer/ConsumerI.h"
#include "../graph/Subgraph.h"


class Problem {
public:
    std::vector<VertexPair> pairs;
    bool solved = false;

public:
    [[nodiscard]] bool is_solved() const {
        return solved;
    }

    [[nodiscard]] auto begin() const {
        return pairs.begin();
    }

    [[nodiscard]] auto end() const {
        return pairs.end();
    }
};

class SelectorI : public ConsumerI {
protected:
    std::shared_ptr<FinderI> finder;
public:
    explicit SelectorI(std::shared_ptr<FinderI> finder_ptr) : finder(std::move(finder_ptr)) {}

    virtual Problem select_problem(Cost k) = 0;

    enum class RecursionType {
        Subgraph,
        VertexPair
    };
    [[nodiscard]] virtual RecursionType recursion_type() const { return RecursionType::Subgraph; }
};

#endif //WEIGHTED_F_FREE_EDGE_EDITING_SELECTORI_H
