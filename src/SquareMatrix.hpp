/**
 * @file SquareMatrix.hpp
 * @author Ryszard Mleczko
 * @date Jan 2025
 * @brief A templated SquareMatrix class for general-purpose mathematical operations.
 */

#include "Matrix.hpp"

namespace rm {
template <typename T>
/**
 * @class SquareMatrix
 * @brief Represents a SquareMatrix derived from Matrix with various methods
 */
class SquareMatrix : public Matrix<T> {
public:

    /**
     * Constructors and destructor
     */

    /**
     * @brief Parametric constructor of SquareMatrix with default value
     * @param rows - number of rows with which SquareMatrix is initialized
     * @param cols - number of columns with which SquareMatrix is initialized
     * @param initValue - default value with which SquareMatrix is initialized
     */
    SquareMatrix(int dim, T initValue = T()) : Matrix<T>(dim, dim, initValue) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension of the SquareMatrix must be greater than zero.");
        }
    }

    /**
     * @brief Parametric constructor of SquareMatrix with initializer list
     * @param rows - number of rows with which SquareMatrix is initialized
     * @param cols - number of columns with which SquareMatrix is initialized
     * @param values - values with which SquareMatrix is initialized
     */
    SquareMatrix(int dim, std::initializer_list<std::initializer_list<T>> values)
        : Matrix<T>(dim, dim, values) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension of the SquareMatrix must be greater than zero.");
        }
    }

    /**
     * @brief Copy constructor of SquareMatrix
     * @param other - SquareMatrix of which is copy being made of 
     */
    SquareMatrix(const SquareMatrix& other) : Matrix<T>(other) {}

    /**
     * @brief Copy constructor of SquareMatrix
     * @param other - Matrix of which is copy being made of 
     */
    SquareMatrix(const Matrix<T>& mat) : Matrix<T>(mat) {
        if (mat.getRows() != mat.getCols()) {
            throw std::invalid_argument("Matrix must be square.");
        }
    }

    /**
     * Matrix specific methods
     */
    
    /**
     * @brief Calculates the cofactor of the SquareMatrix at the specified row and column
     * @param row - the row of the SquareMatrix
     * @param col - the column of the SquareMatrix
     * @return the cofactor of the SquareMatrix at the specified row and column
     */
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

    /**
     * @brief Calculates the determinant of the SquareMatrix
     * @return the determinant of the SquareMatrix
     */
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

    /**
     * @brief Calculates the complement of the SquareMatrix
     * @return the complement of the SquareMatrix
     */
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

    /**
     * @brief Calculates the inverse of the SquareMatrix
     * @return the inverse of the SquareMatrix
     */
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

    /**
     * @brief Raises the SquareMatrix to the power of n
     * @param n - the power to raise the SquareMatrix to
     * @return the SquareMatrix raised to the power of n
     */
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

    /**
     * @brief Calculates the trace of the SquareMatrix
     * @return the trace of the SquareMatrix
     */
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

    /** 
     * @brief Extracts a submatrix from the current matrix
     * @param startRow - the starting row index of the submatrix
     * @param startCol - the starting column index of the submatrix
     * @param numRows - the number of rows in the submatrix
     * @param numCols - the number of columns in the submatrix
     * @return A new Matrix object representing the submatrix
     */
    SquareMatrix subMatrix(int startRow, int startCol, int numRows, int numCols) const {
        if (startRow < 0 || startCol < 0 || startRow + numRows > this->getRows() || startCol + numCols > this->getCols()) {
            throw std::out_of_range("Submatrix dimensions exceed matrix boundaries.");
        }
        
        SquareMatrix result(numRows);
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                result.setVal(i, j, this->getVal(startRow + i, startCol + j));
            }
        }

        return result;
    }

    /** 
     * @brief Transposes the current matrix 
     * @return A new Matrix object representing the transposed matrix
     */
    SquareMatrix transpose() const {
        SquareMatrix<T> transposed(this->getRows());
        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < this->getCols(); ++j) {
                transposed.setVal(j, i, this->getVal(i, j));
            }
        }
        return transposed;
    }

    /**
     * @brief Rotates matrix by 90 degrees clockwise 
     */
    SquareMatrix rotate90() const {
        Matrix result(this->getCols(), this->getRows());

        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < this->getCols(); ++j) {
                result.setVal(j, this->getRows() - 1 - i, this->getVal(i, j));
            }
        }

        return result;
    }

    /**
     * @brief Rotates matrix by 180 degrees clockwise 
     */
    SquareMatrix rotate180() const {
        return this->rotate90().rotate90();
    }

    /**
     * @brief Rotates matrix by 270 degrees clockwise 
     */
    SquareMatrix rotate270() const {
        return this->rotate180().rotate90();
    }

    /**
     * Operators
     */

    /**
     * @brief Division operator for square matrices
     * @param other - the SquareMatrix to divide by
     * @return the result of the division
     */
    SquareMatrix operator/(const SquareMatrix& other) const {
        if (this->getRows() != this->getCols() || other.getRows() != other.getCols() || this->getCols() != other.getRows()) {
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

    /**
     * @brief Division assignment operator for square matrices
     * @param other - the SquareMatrix to divide by
     * @return the result of the division
     * 
     */
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
