#include <iostream>
#include "src/Matrix.hpp"         // Include the Matrix class file

int main() {
    // Create two ComplexNumber objects for testing
    ComplexNumber c1(1.0, 2.0); // 1 + 2i
    ComplexNumber c2(3.0, 4.0); // 3 + 4i

    // 1. Matrix creation using default constructor
    rm::Matrix<ComplexNumber> matrix1(2, 2);
    matrix1.setVal(0, 0, c1);
    matrix1.setVal(0, 1, c2);
    matrix1.setVal(1, 0, ComplexNumber(5.0, 6.0));
    matrix1.setVal(1, 1, ComplexNumber(7.0, 8.0));

    std::cout << "Matrix1 (after setting values):\n" << matrix1 << std::endl;

    // 2. Test Matrix creation using initializer list
    rm::Matrix<ComplexNumber> matrix2(2, 2, {{ComplexNumber(1, 1), ComplexNumber(2, 2)}, {ComplexNumber(3, 3), ComplexNumber(4, 4)}});
    std::cout << "Matrix2 (created with initializer list):\n" << matrix2 << std::endl;

    // 3. Test matrix getVal and setVal methods
    std::cout << "Element at (0, 0) of matrix1: " << matrix1.getVal(0, 0) << std::endl;
    std::cout << "Element at (1, 1) of matrix1: " << matrix1.getVal(1, 1) << std::endl;

    // 4. Test matrix fill method
    matrix2.fill(ComplexNumber(0, 0)); // Fill matrix2 with 0 + 0i
    std::cout << "Matrix2 (after fill with 0):\n" << matrix2 << std::endl;

    // 5. Test randomFill method
    matrix1.randomFill(1, 10); // Fill matrix1 with random complex numbers within range [1, 10]
    std::cout << "Matrix1 (after random fill):\n" << matrix1 << std::endl;

    // 6. Test matrix addition
    rm::Matrix<ComplexNumber> matrixAdd = matrix1 + matrix2;
    std::cout << "Matrix1 + Matrix2:\n" << matrixAdd << std::endl;

    // 7. Test matrix subtraction
    rm::Matrix<ComplexNumber> matrixSub = matrix1 - matrix2;
    std::cout << "Matrix1 - Matrix2:\n" << matrixSub << std::endl;

    // 8. Test matrix multiplication
    rm::Matrix<ComplexNumber> matrixMul = matrix1 * matrix2;
    std::cout << "Matrix1 * Matrix2:\n" << matrixMul << std::endl;

    // 9. Test matrix equality
    rm::Matrix<ComplexNumber> matrixEqual(2, 2);
    matrixEqual.setVal(0, 0, ComplexNumber(1, 1));
    matrixEqual.setVal(0, 1, ComplexNumber(2, 2));
    matrixEqual.setVal(1, 0, ComplexNumber(3, 3));
    matrixEqual.setVal(1, 1, ComplexNumber(4, 4));

    if (matrix2 == matrixEqual) {
        std::cout << "Matrix2 and MatrixEqual are equal.\n";
    } else {
        std::cout << "Matrix2 and MatrixEqual are not equal.\n";
    }

    // 10. Test matrix inequality
    if (matrix1 != matrix2) {
        std::cout << "Matrix1 and Matrix2 are not equal.\n";
    } else {
        std::cout << "Matrix1 and Matrix2 are equal.\n";
    }

    // 11. Test assignment operator
    rm::Matrix<ComplexNumber> matrixAssign = matrix1;
    std::cout << "MatrixAssign (after assignment from Matrix1):\n" << matrixAssign << std::endl;

    // 12. Test matrix multiplication assignment
    matrix1 *= matrix2;
    std::cout << "Matrix1 (after multiplication assignment with Matrix2):\n" << matrix1 << std::endl;

    return 0;
}
