#include "src/DiagonalMatrix.hpp"

int main() {
    try {
        // Matrix initialization with default values
        rm::Matrix<float> m1(3, 3, 0); // 3x3 matrix filled with 0
        std::cout << "Matrix m1 initialized with zeros:\n" << m1;

        // Matrix initialization with specific values
        rm::Matrix<float> m2(3, 3, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
        std::cout << "Matrix m2 initialized with values {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}:\n" << m2;

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

        // Submatrix
        rm::Matrix<float> subM = m1.subMatrix(1, 1, 2, 2);
        std::cout << "Submatrix of m1 (starting at (1,1) with size 2x2):\n" << subM;

        // Add row to matrix
        m1.addRow(1, {1, 2, 3});
        std::cout << "Matrix m1 after adding row at position 2:\n" << m1;

        // Add column to matrix
        m1.addCol(2, {5, 6, 7});
        std::cout << "Matrix m1 after adding column at position 2:\n" << m1;

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
        rm::SquareMatrix<float> sm1(3, 1);  // 3x3 square matrix filled with 1
        std::cout << "SquareMatrix sm1 initialized with 1:\n" << sm1;

        // SquareMatrix initialization with values
        rm::SquareMatrix<float> sm2(3, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
        std::cout << "SquareMatrix sm2 initialized with values {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}:\n" << sm2;

        // SquareMatrix initialization with copy constructor
        rm::SquareMatrix<float> sm3 = sm2;
        std::cout << "SquareMatrix sm3 copied from sm2: \n" << sm3;

        // SquareMatrix determinant
        float det = sm3.det();
        std::cout << "Determinant of sm3: " << det << "\n";

        // SquareMatrix inverse
        rm::SquareMatrix<float> sm4 = sm3.inversed();
        std::cout << "Inverse of sm3:\n" << sm4;

        // SquareMatrix trace
        float trace = sm3.trace();
        std::cout << "Trace of sm3: " << trace << "\n";

        // SquareMatrix power
        rm::SquareMatrix<float> sm5 = sm2;
        sm5.fill(2);
        sm5 = sm5.power(3);
        std::cout << "SquareMatrix sm5 filled with 2 then to the power of 2:\n" << sm5;
        
        // DiagonalMatrix initialization with default values
        rm::DiagonalMatrix<float> dm1(3, 5);  // 3x3 diagonal matrix filled with 5
        std::cout << "DiagonalMatrix dm1 initialized with 5 on diagonal:\n" << dm1;

        // DiagonalMatrix initialization with specific values
        rm::DiagonalMatrix<float> dm2(3, {1.5f, 2.5f, 3.5f});
        std::cout << "DiagonalMatrix dm2 initialized with values {1.5, 2.5, 3.5}:\n" << dm2;

        // DiagonalMatrix copy constructor
        rm::DiagonalMatrix<float> dm3 = dm2;
        std::cout << "DiagonalMatrix dm3 copied from dm2:\n" << dm3;

        // DiagonalMatrix fill with specific value
        dm1.fill(3);
        std::cout << "DiagonalMatrix dm1 after fill(3):\n" << dm1;

        // DiagonalMatrix random fill
        dm1.randomFill(0, 10);
        std::cout << "DiagonalMatrix dm1 after randomFill(0, 10):\n" << dm1;

        // DiagonalMatrix determinant
        float detD = dm1.det();
        std::cout << "Determinant of dm1: " << detD << "\n";

        // DiagonalMatrix inverse
        rm::DiagonalMatrix<float> dm4 = dm1.inversed();
        std::cout << "Inverse of dm1:\n" << dm4;

        // DiagonalMatrix trace
        float traceD = dm1.trace();
        std::cout << "Trace of dm1: " << traceD << "\n";

        // DiagonalMatrix isIdentity
        if (dm1.isIdentity()) {
            std::cout << "dm1 is an identity matrix\n";
        } else {
            std::cout << "dm1 is not an identity matrix\n";
        }

        // DiagonalMatrix getDiagonal
        float* diagonal = dm1.getDiagonal();
        std::cout << "Diagonal of dm1: ";
        for (int i = 0; i < dm1.getRows(); ++i) {
            std::cout << diagonal[i] << " ";
        }
        std::cout << "\n";
        
        // DiagonalMatrix setDiagonal
        dm1.setDiagonal({1, 2, 3});
        std::cout << "DiagonalMatrix dm1 after setDiagonal({1, 2, 3}):\n" << dm1;

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    return 0;
}
