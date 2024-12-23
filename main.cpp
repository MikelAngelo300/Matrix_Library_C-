#include "src/Matrix.hpp"

int main() {
    // Making 3x2 Matrix 
    Matrix<int> mat(3, 2, 5);
    std::cout << "Macierz 2x2:\n" << mat;
    
    // Making 3x3 Matrix
    Matrix<double> mat2(3, 3, 1.1);
    std::cout << "Macierz 3x3:\n" << mat2;

    return 0;
}
