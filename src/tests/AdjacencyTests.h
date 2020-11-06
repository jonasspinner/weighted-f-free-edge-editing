#ifndef WEIGHTED_F_FREE_EDGE_EDITING_ADJACENCYTESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_ADJACENCYTESTS_H


#include "test_utils.h"

class AdjacencyTests {

public:
    static void empty_row() {
        Graph::AdjRow row(100);
        Graph::RowVertices row_range(row);

        std::vector<Vertex> actual(row_range.begin(), row_range.end());
        expect("empty row", {}, actual);
    }

    static void first_block_index() {
        Vertex pos = Graph::AdjRow::bits_per_block;
        Graph::AdjRow row(pos + 10);
        row[pos] = true;
        Graph::RowVertices row_range(row);

        std::vector<Vertex> actual(row_range.begin(), row_range.end());
        expect("first block index", {pos}, actual);
    }

    static void last_block_index() {
        Vertex pos = Graph::AdjRow::bits_per_block - 1;
        Graph::AdjRow row(pos + 10);
        row[pos] = true;
        Graph::RowVertices row_range(row);

        std::vector<Vertex> actual(row_range.begin(), row_range.end());
        expect("last block index", {pos}, actual);
    }

    void run() {
        empty_row();
        first_block_index();
        last_block_index();
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_ADJACENCYTESTS_H
