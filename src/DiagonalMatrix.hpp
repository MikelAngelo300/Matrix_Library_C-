#include "SquareMatrix.hpp"

namespace rm {
template <typename T>
requires std::is_arithmetic_v<T>
class DiagonalMatrix : public SquareMatrix<T> {
public:
    
    DiagonalMatrix(int dim, T initValue = T()) : SquareMatrix<T>(dim, T()) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension of the matrix must be greater than zero.");
        }

        for (int i=0; i<dim;i++) {
            this->m_data[i][i]=initValue;
        }
    }

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

    DiagonalMatrix(const SquareMatrix<T>& other) : SquareMatrix<T>(other) {
        for (int i = 0; i < this->m_rows; ++i) {
            for (int j = 0; j < this->m_cols; ++j) {
                if (i != j && this->m_data[i][j] != T()) {
                    throw std::invalid_argument("Input matrix is not diagonal.");
                }
            }
        }
    }

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

    void setDiagonal(std::initializer_list<T> values) {
        if (values.size() != static_cast<size_t>(this->m_rows)) {
            throw std::invalid_argument("Number of diagonal elements must match the dimension of the matrix.");
        }

        auto it = values.begin();
        for (int i = 0; i < this->m_rows; ++i) {
            this->m_data[i][i] = *it++;
        }
    }
    
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