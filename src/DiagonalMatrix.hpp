/**
 * @file DiagonalMatrix.hpp
 * @brief Implementation of a templated DiagonalMatrix class derived from SquareMatrix.
 *
 * This class represents diagonal matrices, where only the diagonal elements are non-zero.
 * It builds upon the SquareMatrix class, optimizing for diagonal-specific operations and
 * ensuring constraints are maintained during initialization and modification.
 * 
 * Features include:
 * - Constructors to initialize with a uniform value or a list of diagonal elements.
 * - Validation of diagonal matrix properties for input matrices.
 * - Extraction and modification of diagonal elements (`getDiagonal`, `setDiagonal`).
 * - Identity matrix check (`isIdentity`).
 * - Runtime validation for matrix dimensions and diagonal-specific constraints.
 * 
 * The DiagonalMatrix class is ideal for cases where mathematical operations on diagonal matrices
 * are required, leveraging the efficiency of operations limited to diagonal elements.
 */

#include "SquareMatrix.hpp"

namespace rm {
template <typename T>
requires std::is_arithmetic_v<T>
class DiagonalMatrix : public SquareMatrix<T> {
public:
    // Constructor for DiagonalMatrix with default value
    DiagonalMatrix(int dim, T initValue = T()) : SquareMatrix<T>(dim, T()) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension of the matrix must be greater than zero.");
        }

        for (int i = 0; i < dim; i++) {
            this->setVal(i,i, initValue);
        }
    }

    // Constructor for DiagonalMatrix with initializer list
    DiagonalMatrix(int dim, std::initializer_list<T> values)
        : SquareMatrix<T>(dim) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension of the matrix must be greater than zero.");
        }
        if (values.size() != static_cast<size_t>(dim)) {
            throw std::invalid_argument("Number of diagonal elements must match the dimension of the matrix.");
        }

        auto it = values.begin();
        for (int i = 0; i < dim; ++i) {
            this->setVal(i,i, *it++);
        }
    }

    // Copy constructor
    DiagonalMatrix(const DiagonalMatrix<T>& other) : SquareMatrix<T>(other) {
        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < this->getCols(); ++j) {
                if (i != j && this->getVal(i,j) != T()) {
                    throw std::invalid_argument("Input matrix is not diagonal.");
                }
            }
        }
    }

    // Getter for diagonal
    T* getDiagonal() const {
    if (this->getRows() != this->getCols()) {
        throw std::logic_error("Matrix is not square, cannot extract diagonal.");
    }

    T* diagonal = new T[this->getRows()];
    for (int i = 0; i < this->getRows(); ++i) {
        diagonal[i] = this->getVal(i, i);  // Extract diagonal elements
    }

    return diagonal;
}

    // Setter for diagonal
    void setDiagonal(std::initializer_list<T> values) {
        if (values.size() != static_cast<size_t>(this->getRows())) {
            throw std::invalid_argument("Number of diagonal elements must match the dimension of the matrix.");
        }

        auto it = values.begin();
        for (int i = 0; i < this->getRows(); ++i) {
            this->setVal(i,i, *it++);
        }
    }

    // Checks whether the matrix is an identity matrix
    bool isIdentity() const {
        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < this->getCols(); ++j) {
                if (i != j && this->getVal(i,j) != T()) {
                    return false;
                }
                if (i == j && this->getVal(i,j) != T(1)) {
                    return false;
                }
            }
        }
        return true;
    }

    void fill(T value) override {
        for (int i = 0; i < this->getRows(); ++i) {
            this->setVal(i,i, value);
        }
    }

    // Fill the diagonal with random values within a specified range
    void randomFill(int min, int max) override {
        if (min > max) {
            throw std::invalid_argument("min cannot be greater than max.");
        }

        static bool seeded = false;
        if (!seeded) {
            srand(static_cast<unsigned int>(time(0)));
            seeded = true;
        }

        for (int i = 0; i < this->getRows(); ++i) {
            this->setVal(i,i, rand() % (max - min + 1) + min);
        }
    }

    // Fill the diagonal with Gaussian distribution of random values
    void randomFillGaussian(T mean, T stddev) override {
        static_assert(std::is_floating_point<T>::value, "T must be a floating-point type.");
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<T> dist(mean, stddev);

        for (int i = 0; i < this->getRows(); ++i) {
            this->setVal(i,i,dist(gen));
        }
    }

    // Fill the matrix with 0
    void zeros() override {
        this->fill(0);
    }

    // Fill the diagonal with 1
    void ones() override {
        this->fill(1);
    }
    
};
}
