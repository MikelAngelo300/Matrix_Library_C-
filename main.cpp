#include "src/Matrix.hpp"

int main() {
    // Tworzymy macierz 2x2 typu int
    Matrix<int> mat(3, 2, 5);
    std::cout << "Macierz 2x2:\n" << mat;
    
    // Tworzymy macierz 3x3 typu double
    Matrix<double> mat2(3, 3, 1.1);
    std::cout << "Macierz 3x3:\n" << mat2;

    return 0;
}
