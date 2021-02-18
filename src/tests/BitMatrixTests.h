
#ifndef WEIGHTED_F_FREE_EDGE_EDITING_BITMATRIXTESTS_H
#define WEIGHTED_F_FREE_EDGE_EDITING_BITMATRIXTESTS_H

#include "../graph/BitMatrix.h"

class BitMatrixTests {
public:

    void test(const dynamic_bitset::BitMatrix<> &matrix) {
        auto row = matrix.row(2);
        std::cout << row[3];
    }

    void run() {
        dynamic_bitset::BitMatrix<> matrix(10);
        dynamic_bitset::BitMatrix<> matrix2(10);

        auto row2 = matrix2.mut_row(2);
        row2[1] = true;
        matrix.mut_row(2)[3] = true;
        matrix.mut_row(2)[4] = true;
        matrix[{3, 5}] = true;

        matrix.mut_row(3)[1] = true;
        matrix.mut_row(4).flip(2);
        matrix.mut_row(5).flip();
        matrix.mut_row(6).reset(4);
        matrix.mut_row(7).reset();
        matrix.mut_row(8).set(6);
        matrix.mut_row(9).set();

        matrix.mut_row(9) -= matrix.row(1);
        matrix.mut_row(8) ^= matrix.row(2);
        matrix.mut_row(7) |= matrix.row(3);
        matrix.mut_row(6) &= matrix.row(4);

        matrix.mut_row(5).flip();
        matrix.mut_row(5) -= matrix.row(2);
        matrix2.mut_row(2) = matrix.row(8);

        std::cout << matrix << "\n";
        std::cout << matrix2 << "\n";

        dynamic_bitset::BitRow<> row(10);
        row = matrix.row(2);

        std::cout << row << "\n";

        std::cout << std::vector<Vertex>(matrix.row(9).indices().begin(), matrix.row(9).indices().end()) << "\n";
        std::cout << std::vector<Vertex>(row.indices().begin(), row.indices().end()) << "\n";


        dynamic_bitset::BitRow<> long_row(10000);
        long_row.set(8129);
        long_row.set(4359);
        long_row.set(1583);

        std::cout << std::vector<Vertex>(long_row.indices().begin(), long_row.indices().end()) << "\n";
    }
};


#endif //WEIGHTED_F_FREE_EDGE_EDITING_BITMATRIXTESTS_H
