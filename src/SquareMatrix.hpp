/**
 * @file SquareMatrix.hpp
 * @author Ryszard Mleczko
 * @date Jan 2025
 * @brief A templated square matrix class for general-purpose mathematical operations.
 * @tparam T The data type of the matrix elements (must be arithmetic or satisfy custom requirements).
 */

#include "Matrix.hpp"

namespace rm {
template <typename T>
class SquareMatrix : public Matrix<T> {
public:
    // Constructor for SquareMatrix with defualt value
    SquareMatrix(int dim, T initValue = T()) : Matrix<T>(dim, dim, initValue) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension of the matrix must be greater than zero.");
        }
    }

    // Constructor for SquareMatrix with initializer list
    SquareMatrix(int dim, std::initializer_list<std::initializer_list<T>> values)
        : Matrix<T>(dim, dim, values) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension of the matrix must be greater than zero.");
        }
    }

    // Copy constructor
    SquareMatrix(const SquareMatrix& other) : Matrix<T>(other) {}

    // Constructor to initialize from Matrix<T>
    SquareMatrix(const Matrix<T>& mat) : Matrix<T>(mat) {
        if (mat.getRows() != mat.getCols()) {
            throw std::invalid_argument("Matrix must be square.");
        }
    }


    // Function to calculate the cofactor
    T cofactor(int row, int col) const {
        if (this->getRows() != this->getCols()) {
            throw std::invalid_argument("Matrix must be square to compute its cofactor.");
        }

        SquareMatrix<T> minorMatrix = this->subMatrix(0, 0, this->getRows(), this->getCols());
        minorMatrix.removeRow(row);
        minorMatrix.removeColumn(col);

        T determinant = minorMatrix.det();

        return ((row + col) % 2 == 0 ? 1 : -1) * determinant;
    }

    // Function to calculate the determinant
    T det() const {
        if (this->getRows() != this->getCols()) {
            throw std::invalid_argument("Matrix must be square to compute its determinant.");
        }

        if (this->getRows() == 2) {
            return this->getVal(0, 0) * this->getVal(1, 1) - this->getVal(0, 1) * this->getVal(1, 0);
        }

        T determinant = 0;
        for (int i = 0; i < this->getCols(); ++i) {
            determinant += (i % 2 == 0 ? 1 : -1) * this->getVal(0, i) * this->cofactor(0, i);
        }
        return determinant;
    }

    // Function to calculate the cofactor matrix
    SquareMatrix complement() const {
        if (this->det() == 0) {
            throw std::invalid_argument("Matrix is singular (determinant is 0)!");
        }

        SquareMatrix result(this->getRows());

        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < this->getCols(); ++j) {
                result.setVal(i, j, cofactor(i, j) * ((i + j) % 2 == 0 ? 1 : -1));
            }
        }

        return result;
    }

    // Function to calculate the inverse of the matrix
    SquareMatrix<T> inversed() const {
        if (this->det() == 0) {
            throw std::invalid_argument("Matrix is singular and cannot be inverted.");
        }

        Matrix<T> comp = this->complement();

        Matrix<T> trans = comp.transpose();

        if (trans.getRows() != this->getRows() || trans.getCols() != this->getCols()) {
            throw std::invalid_argument("Transpose dimensions are inconsistent.");
        }

        SquareMatrix<T> squareTrans(trans.getRows(), trans.getCols());
        for (int i = 0; i < trans.getRows(); ++i) {
            for (int j = 0; j < trans.getCols(); ++j) {
                squareTrans.setVal(i, j, trans.getVal(i, j));
            }
        }

        T determinant = this->det();
        if (determinant != 0) {
            squareTrans.scale(1 / determinant);
        } else {
            throw std::invalid_argument("Matrix is singular, cannot scale.");
        }

        return squareTrans;
    }

    // Function to calculate the power of the matrix
    SquareMatrix<T> power(int n) const {
        if (n < 0) {
            throw std::invalid_argument("Exponent must be non-negative.");
        }

        SquareMatrix<T> result(this->getRows(), T(0));
        for (int i = 0; i < this->getRows(); ++i) {
            result.setVal(i, i, T(1)); 
        }

        SquareMatrix<T> base = *this;

        for (int i = 0; i < n; ++i) {
            result *= base; 
        }

        return result;
    }

    // Function to calculate the trace of the matrix
    T trace() const {
        if (this->getRows() != this->getCols()) {
            throw std::invalid_argument("Matrix must be square to calculate the trace.");
        }

        T result = 0;
        for (int i = 0; i < this->getRows(); ++i) {
            result += this->getVal(i, i);
        }

        return result;
    }

    // Square matrix division
    SquareMatrix operator/(const SquareMatrix& other) const {
        if (this->getRows() != this->getCols() || other.m_rows != other.m_cols || this->getCols() != other.m_rows) {
            throw std::invalid_argument("Matrices must be square and of the same size for division.");
        }
        SquareMatrix inverse = other.inversed();
        SquareMatrix result(this->getRows(), this->getCols());

        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < this->getCols(); ++j) {
                T sum = 0;
                for (int k = 0; k < this->getCols(); ++k) {
                    sum += this->m_data[i][k] * inverse.m_data[k][j];
                }
                result.m_data[i][j] = sum;
            }
        }

        return result;
    }

    // Square matrix division assignment
    SquareMatrix& operator/=(const SquareMatrix& other) {
        if (this->getRows() != this->getCols() || other.m_rows != other.m_cols || this->getCols() != other.m_rows) {
            throw std::invalid_argument("Matrices must be square and of the same size for division.");
        }
        SquareMatrix inverse = other.inversed();
        SquareMatrix result(this->getRows(), this->getCols());

        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < this->getCols(); ++j) {
                T sum = 0;
                for (int k = 0; k < this->getCols(); ++k) {
                    sum += this->m_data[i][k] * inverse.m_data[k][j];
                }
                result.m_data[i][j] = sum;
            }
        }

        *this = std::move(result);
        return *this;
    }
};
}
