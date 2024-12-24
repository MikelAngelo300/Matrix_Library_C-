#include "src/Matrix.hpp"
// Function to test Matrix class with complex numbers
int main() {
    try {
        // Create a 2x2 matrix of ComplexNumber objects
        rm::Matrix<ComplexNumber> mat1(2, 2, ComplexNumber(1, 2));  // Constructor with default value
        std::cout << "Matrix 1 (initialized with Complex(1,2)):\n" << mat1 << std::endl;

        // Use randomFill to fill the matrix with random values (real numbers here)
        mat1.randomFill(1, 5);
        std::cout << "Matrix 1 after random fill:\n" << mat1 << std::endl;

        // Create another matrix with complex numbers directly
        rm::Matrix<ComplexNumber> mat2(2, 2, { { ComplexNumber(1, 1), ComplexNumber(2, 2) },
                                               { ComplexNumber(3, 3), ComplexNumber(4, 4) } });
        std::cout << "Matrix 2 (initialized with Complex numbers):\n" << mat2 << std::endl;

        // Matrix addition
        auto matAdd = mat1 + mat2;
        std::cout << "Matrix 1 + Matrix 2:\n" << matAdd << std::endl;

        // Matrix subtraction
        auto matSub = mat1 - mat2;
        std::cout << "Matrix 1 - Matrix 2:\n" << matSub << std::endl;

        // Matrix multiplication
        auto matMul = mat1 * mat2;
        std::cout << "Matrix 1 * Matrix 2:\n" << matMul << std::endl;

        // Matrix addition assignment
        mat1 += mat2;
        std::cout << "Matrix 1 after Matrix 1 += Matrix 2:\n" << mat1 << std::endl;

        // Matrix subtraction assignment
        mat1 -= mat2;
        std::cout << "Matrix 1 after Matrix 1 -= Matrix 2:\n" << mat1 << std::endl;

        // Matrix multiplication assignment
        mat1 *= mat2;
        std::cout << "Matrix 1 after Matrix 1 *= Matrix 2:\n" << mat1 << std::endl;

        // Test the equality operator (it should not be equal after multiplication)
        std::cout << "Matrix 1 == Matrix 2: " << (mat1 == mat2 ? "True" : "False") << std::endl;

        // Test the inequality operator
        std::cout << "Matrix 1 != Matrix 2: " << (mat1 != mat2 ? "True" : "False") << std::endl;

        // Test setting and getting individual values
        mat1.setVal(0, 0, ComplexNumber(10, 10));
        std::cout << "Matrix 1 after setting (0, 0) to Complex(10,10):\n" << mat1 << std::endl;
        ComplexNumber value = mat1.getVal(0, 0);
        std::cout << "Matrix 1[0][0] = " << value << std::endl;

        // Test the copy constructor
        rm::Matrix<ComplexNumber> mat3(mat2); // Copy constructor
        std::cout << "Matrix 3 (copy of Matrix 2):\n" << mat3 << std::endl;

        // Test the assignment operator
        rm::Matrix<ComplexNumber> mat4(2, 2);
        mat4 = mat1; // Assignment operator
        std::cout << "Matrix 4 (assigned from Matrix 1):\n" << mat4 << std::endl;

        // Test destructor by creating a large matrix
        {
            rm::Matrix<int> largeMat(100, 100);
            largeMat.randomFill(1, 10);
            std::cout << "Large matrix created and filled." << std::endl;
        } // Destructor will be called here when going out of scope

    } catch (const std::exception& e) {
        std::cout << "Caught exception: " << e.what() << std::endl;
    }

    return 0;
}
