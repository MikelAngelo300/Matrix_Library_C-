/**
 * @file SquareMatrix.hpp
 * @brief Implementation of a templated SquareMatrix class derived from the Matrix base class.
 *
 * This class provides operations specific to square matrices, including determinant calculation,
 * inversion, matrix power, trace, and more. It is designed to handle mathematical operations
 * commonly used in linear algebra while ensuring robust error checking for invalid operations.
 * 
 * Features include:
 * - Determinant computation (`det`).
 * - Matrix inversion (`inversed`).
 * - Cofactor and cofactor matrix calculation (`cofactor`, `complement`).
 * - Matrix power (`power`).
 * - Trace calculation (`trace`).
 * - Square matrix division (`operator/` and `operator/=`).
 * - Constructors for initialization with dimensions, values, or another matrix.
 * 
 * This class ensures dimension validation and provides runtime exceptions for invalid operations.
 */

#include "Matrix.hpp"

namespace rm {
template <typename T>
class SquareMatrix : public Matrix<T> {
public:
    // Constructor with size validation for square matrix
    SquareMatrix(int dim, T initValue = T()) : Matrix<T>(dim, dim, initValue) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension of the matrix must be greater than zero.");
        }
    }

    // Constructor with initializer list
    SquareMatrix(int dim, std::initializer_list<std::initializer_list<T>> values)
        : Matrix<T>(dim, dim, values) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension of the matrix must be greater than zero.");
        }
    }

    // Copy constructor
    SquareMatrix(const SquareMatrix& other) : Matrix<T>(other) {}

    // Function to calculate the cofactor
    T cofactor(int row, int col) const {
        if (this->m_rows != this->m_cols) {
            throw std::invalid_argument("Matrix must be square to compute its cofactor.");
        }

        // Create a submatrix excluding the specified row and column
        SquareMatrix minorMatrix = this->subMatrix(0, 0, this->m_rows, this->m_cols).removeRow(row).removeColumn(col);

        // Compute the determinant of the minor matrix
        T determinant = minorMatrix.det();

        // Adjust the sign based on the cofactor position
        return ((row + col) % 2 == 0 ? 1 : -1) * determinant;
    }


    // Function to calculate the determinant
    T det() const {
        if (this->m_rows != this->m_cols) {
            throw std::invalid_argument("Matrix must be square to compute its determinant.");
        }

        if (this->m_rows == 2) {
            return this->getVal(0, 0) * this->getVal(1, 1) - this->getVal(0, 1) * this->getVal(1, 0);
        }

        T determinant = 0;
        for (int i = 0; i < this->m_cols; ++i) {
            determinant += (i % 2 == 0 ? 1 : -1) * this->getVal(0, i) * cofactor(0, i);
        }
        return determinant;
    }

    // Function to calculate the cofactor matrix
    SquareMatrix complement() const {
        if (this->det() == 0) {
            throw std::invalid_argument("Matrix is singular (determinant is 0)!");
        }

        SquareMatrix result(this->m_rows);

        for (int i = 0; i < this->m_rows; ++i) {
            for (int j = 0; j < this->m_cols; ++j) {
                // Set cofactor with the proper sign
                result.setVal(i, j, cofactor(i, j) * ((i + j) % 2 == 0 ? 1 : -1));
            }
        }

        return result;
    }

    // Function to calculate the inverse of the matrix
    SquareMatrix inversed() const {
        if (this->det() == 0) {
            throw std::invalid_argument("Matrix is singular and cannot be inverted.");
        }

        SquareMatrix comp = this->complement();
        SquareMatrix trans = comp.transpose();  // Transpose the cofactor matrix
        T determinant = this->det();

        // Use the inherited scale function to scale the transposed matrix
        SquareMatrix result = trans.scale(1 / determinant);  // Scale by 1/determinant

        return result;
    }

    // Function to calculate the power of the matrix
    SquareMatrix power(int n) const {
        if (n < 0) {
            throw std::invalid_argument("Exponent must be non-negative.");
        }
        
        // Initialize the result as the identity matrix
        SquareMatrix result(this->m_rows, T(1));  // Identity matrix of size m_rows
        SquareMatrix base = *this;

        // Multiply result by base (matrix) n times
        for (int i = 0; i < n; ++i) {
            result *= base;  // Multiply result by base (this) in each iteration
        }

        return result;
    }

    // Function to calculate the trace of the matrix
    T trace() const {
        if (this->m_rows != this->m_cols) {
            throw std::invalid_argument("Matrix must be square to calculate the trace.");
        }

        T result = 0;
        for (int i = 0; i < this->m_rows; ++i) {
            result += this->getVal(i, i);  // Only sum diagonal elements
        }

        return result;
    }

    // Square matrix division
    SquareMatrix operator/(const SquareMatrix& other) const {
        if (m_rows != m_cols || other.m_rows != other.m_cols || m_cols != other.m_rows) {
            throw std::invalid_argument("Matrices must be square and of the same size for division.");
        }
        SquareMatrix inverse = other.inversed();
        SquareMatrix result(m_rows, m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                T sum = 0;
                for (int k = 0; k < m_cols; ++k) {
                    sum += m_data[i][k] * inverse.m_data[k][j];
                }
                result.m_data[i][j] = sum;
            }
        }

        return result;
    }

    // Square matrix division assignment
    SquareMatrix& operator/=(const SquareMatrix& other) {
        if (m_rows != m_cols || other.m_rows != other.m_cols || m_cols != other.m_rows) {
            throw std::invalid_argument("Matrices must be square and of the same size for division.");
        }
        SquareMatrix inverse = other.inversed();
        SquareMatrix result(m_rows, m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                T sum = 0;
                for (int k = 0; k < m_cols; ++k) {
                    sum += m_data[i][k] * inverse.m_data[k][j];
                }
                result.m_data[i][j] = sum;
            }
        }

        *this = std::move(result);
        return *this;
    }


    

};
}
