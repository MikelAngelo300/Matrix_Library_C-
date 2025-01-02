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
    // Constructor for DiagonalMatrix with defualt value
    DiagonalMatrix(int dim, T initValue = T()) : SquareMatrix<T>(dim, T()) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension of the matrix must be greater than zero.");
        }

        for (int i=0; i<dim;i++) {
            this->m_data[i][i]=initValue;
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
            this->m_data[i][i] = *it++;
        }
    }
    // Copy constructor
    DiagonalMatrix(const DiagonalMatrix<T>& other) : SquareMatrix<T>(other) {
        for (int i = 0; i < this->m_rows; ++i) {
            for (int j = 0; j < this->m_cols; ++j) {
                if (i != j && this->m_data[i][j] != T()) {
                    throw std::invalid_argument("Input matrix is not diagonal.");
                }
            }
        }
    }


    // getter for diagonal
    std::initializer_list<T> getDiagonal() const {
        if (this->m_rows != this->m_cols) {
            throw std::logic_error("Matrix is not square, cannot extract diagonal.");
        }

        T diagonal[this->m_rows]; 
        for (int i = 0; i < this->m_rows; ++i) {
            diagonal[i] = this->m_data[i][i]; 
        }

        return {diagonal, diagonal + this->m_rows}; 
    }

    // Setter for diagonal
    void setDiagonal(std::initializer_list<T> values) {
        if (values.size() != static_cast<size_t>(this->m_rows)) {
            throw std::invalid_argument("Number of diagonal elements must match the dimension of the matrix.");
        }

        auto it = values.begin();
        for (int i = 0; i < this->m_rows; ++i) {
            this->m_data[i][i] = *it++;
        }
    }
    
    // Checks whether there are only ones in the diagonal
    bool isIdentity() const {
    for (int i = 0; i < this->m_rows; ++i) {
        for (int j = 0; j < this->m_cols; ++j) {
            if (i != j && this->m_data[i][j] != T()) {
                return false;
            }
            if (i == j && this->m_data[i][j] != T(1)) {
                return false;
            }
        }
    }
    return true;
}

};
}