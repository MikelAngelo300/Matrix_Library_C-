#include "src/DiagonalMatrix.hpp"

int main() {
    try {
        // Matrix initialization with default values
        rm::Matrix<float> m1(3, 3, 0); // 3x3 matrix filled with 0
        std::cout << "Matrix m1 initialized with zeros:\n" << m1;

        // Matrix initialization with specific values
        rm::Matrix<float> m2(2, 2, {{1, 2}, {3, 4}});
        std::cout << "Matrix m2 initialized with values {{1, 2}, {3, 4}}:\n" << m2;

        // Matrix copy constructor
        rm::Matrix<float> m3 = m2;
        std::cout << "Matrix m3 copied from m2:\n" << m3;

        // Matrix fill with specific value
        m1.fill(5);
        std::cout << "Matrix m1 after fill(5):\n" << m1;

        // Matrix random fill
        m1.randomFill(0, 10);
        std::cout << "Matrix m1 after randomFill(1, 10):\n" << m1;

        
        // Matrix sum and mean
        std::cout << "Matrix m1 sum: " << m1.sum() << "\n";
        std::cout << "Matrix m1 mean: " << m1.mean() << "\n";

        // Matrix min and max
        std::cout << "Matrix m1 min: " << m1.min() << "\n";
        std::cout << "Matrix m1 max: " << m1.max() << "\n";

        // Matrix normalize
        m1.normalize();
        std::cout << "Matrix m1 after normalization:\n" << m1;

        // Matrix transpose
        rm::Matrix<float> m4 = m1.transpose();
        std::cout << "Matrix m4 (transpose of m1):\n" << m4;

        // Matrix rotation (90, 180, 270 degrees)
        rm::Matrix<float> m5 = m1.rotate90();
        std::cout << "Matrix m5 after rotate90:\n" << m5;
        rm::Matrix<float> m6 = m1.rotate180();
        std::cout << "Matrix m6 after rotate180:\n" << m6;
        rm::Matrix<float> m7 = m1.rotate270();
        std::cout << "Matrix m7 after rotate270:\n" << m7;

        // Matrix addition
        rm::Matrix<float> m8 = m1 + m2;
        std::cout << "Matrix m8 after m1 + m2:\n" << m8;

        // Matrix subtraction
        rm::Matrix<float> m9 = m1 - m2;
        std::cout << "Matrix m9 after m1 - m2:\n" << m9;

        // Matrix multiplication
        rm::Matrix<float> m10 = m1 * m2;
        std::cout << "Matrix m10 after m1 * m2:\n" << m10;

        // Matrix scaling
        m1.scale(2);
        std::cout << "Matrix m1 after scaling by 2:\n" << m1;

        // Matrix resizing
        m1.resize(4, 4, 1);
        std::cout << "Matrix m1 after resizing to 4x4:\n" << m1;

        // Submatrix
        rm::Matrix<float> subM = m1.subMatrix(1, 1, 2, 2);
        std::cout << "Submatrix of m1 (starting at (1,1) with size 2x2):\n" << subM;

        // Add row to matrix
        m1.addRow(2, {1, 2, 3, 4});
        std::cout << "Matrix m1 after adding row at position 2:\n" << m1;

        // Add column to matrix
        m1.addCol(3, {5, 6, 7, 8});
        std::cout << "Matrix m1 after adding column at position 3:\n" << m1;

        // Remove row
        m1.removeRow(1);
        std::cout << "Matrix m1 after removing row 1:\n" << m1;

        // Remove column
        m1.removeColumn(2);
        std::cout << "Matrix m1 after removing column 2:\n" << m1;

        // Equality operator
        if (m1 == m2) {
            std::cout << "m1 is equal to m2\n";
        } else {
            std::cout << "m1 is not equal to m2\n";
        }

        // Inequality operator
        if (m1 != m2) {
            std::cout << "m1 is not equal to m2\n";
        }

        // SquareMatrix initialization
        rm::SquareMatrix<float> sm1(3, 1);
        std::cout << "SquareMatrix sm1 initialized with 1:\n" << sm1;

        // SquareMatrix initialization with values
        rm::SquareMatrix<float> sm2(3, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
        std::cout << "SquareMatrix sm2 initialized with values {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}:\n" << sm2;

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    return 0;
}
