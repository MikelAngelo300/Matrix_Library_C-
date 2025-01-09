/**
 * @file Matrix.hpp
 * @author Ryszard Mleczko
 * @date Jan 2025
 * @brief A templated DiagonalMatrix class for general-purpose mathematical operations.
 */


#include "SquareMatrix.hpp"

namespace rm {
template <typename T>
requires std::is_arithmetic_v<T>
class DiagonalMatrix : public SquareMatrix<T> {
public:
    
    /**
     * Constructors and destructor
     */

    /**
     * @brief Parametric constructor of DiagonalMatrix with default value
     * @param dim - number of rows with which DiagonalMatrix is initialized
     * @param initValue - default value with which DiagonalMatrix is initialized
     */
    DiagonalMatrix(int dim, T initValue = T()) : SquareMatrix<T>(dim, T()) {
        if (dim <= 0) {
            throw std::invalid_argument("Dimension of the matrix must be greater than zero.");
        }

        for (int i = 0; i < dim; i++) {
            this->setVal(i,i, initValue);
        }
    }

    /**
     * @brief Parametric constructor of DiagonalMatrix with initializer list
     * @param dim - number of rows with which DiagonalMatrix is initialized
     * @param values - values with which DiagonalMatrix is initialized
     */
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

    /**
     * @brief Copy constructor of DiagonalMatrix
     * @param other - DiagonalMatrix of which is copy being made of
     */
    DiagonalMatrix(const DiagonalMatrix<T>& other) : SquareMatrix<T>(other) {
        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < this->getCols(); ++j) {
                if (i != j && this->getVal(i,j) != T()) {
                    throw std::invalid_argument("Input matrix is not diagonal.");
                }
            }
        }
    }

    /**
     * @brief Copy constructor of DiagonalMatrix from SquareMatrix
     * @param other - SquareMatrix of which is copy being made of
     */
    DiagonalMatrix(const SquareMatrix<T>& other) : SquareMatrix<T>(other) {
        if (this->getRows() != this->getCols()) {
            throw std::invalid_argument("Matrix must be square to be converted to diagonal.");
        }

        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < this->getCols(); ++j) {
                if (i != j && this->getVal(i,j) != T()) {
                    throw std::invalid_argument("Input matrix is not diagonal.");
                }
            }
        
        }
    }

    /**
     * @brief Copy constructor of DiagonalMatrix from Matrix
     * @param other - Matrix of which is copy being made of
     */
    DiagonalMatrix(const Matrix<T>& other) : SquareMatrix<T>(other) {
        if (this->getRows() != this->getCols()) {
            throw std::invalid_argument("Matrix must be square to be converted to diagonal.");
        }

        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < this->getCols(); ++j) {
                if (i != j && this->getVal(i,j) != T()) {
                    throw std::invalid_argument("Input matrix is not diagonal.");
                }
            }
        }
    }

    /**
     * Getters
     */

    /**
     * @brief Gets the diagonal of the matrix
     * @return the diagonal of the matrix
     */
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

    /**
     * Setters
     */

    /**
     * @brief Sets the diagonal of the matrix
     * @param values - values to set the diagonal of the matrix
     */
    void setDiagonal(std::initializer_list<T> values) {
        if (values.size() != static_cast<size_t>(this->getRows())) {
            throw std::invalid_argument("Number of diagonal elements must match the dimension of the matrix.");
        }

        auto it = values.begin();
        for (int i = 0; i < this->getRows(); ++i) {
            this->setVal(i,i, *it++);
        }
    }

    /**
     * @brief Sets value at specific location
     * @param i - row number 
     * @param j - column number
     * @param newVal - new Value that is being set 
     */
    void setVal(int i, int j, T newVal) override {
        if (i != j) {
            throw std::invalid_argument("Only diagonal elements can be set in a DiagonalMatrix.");
        }
        SquareMatrix<T>::setVal(i, j, newVal);
    }

    /**
     * Operations
     */

    /**
     * @brief Checks if the DiagonalMatrix is an identity matrix
     * @return the boolean value indicating if the DiagonalMatrix is an identity matrix
     */
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

    /**
     * Fillers
     */

    /**
     * @brief Fills the diagonal with a specific value
     * @param value - the value to fill the diagonal with
     */
    void fill(T value) override {
        for (int i = 0; i < this->getRows(); ++i) {
            this->setVal(i,i, value);
        }
    }

    /**
     * @brief Fills the diagonal with random values within a specified range
     * @param min - minimal range value
     * @param max - maximal range value
     */
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

    /**
     * @brief Fills the diagonal with random values from a Gaussian distribution
     * @param mean - mean of the Gaussian distribution
     * @param stddev - standard deviation of the Gaussian distribution
     */
    void randomFillGaussian(T mean, T stddev) override {
        static_assert(std::is_floating_point<T>::value, "T must be a floating-point type.");
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<T> dist(mean, stddev);

        for (int i = 0; i < this->getRows(); ++i) {
            this->setVal(i,i,dist(gen));
        }
    }

    /**
     * @brief Fills the diagonal with zeros
     */
    void zeros() override {
        this->fill(0);
    }

    /**
     * @brief Fills the diagonal with ones
     */
    void ones() override {
        this->fill(1);
    }
    
    /**
     * Matrix specific methods
     */

    /** 
     * @brief Extracts a submatrix from the current matrix
     * @param startRow - the starting row index of the submatrix
     * @param startCol - the starting column index of the submatrix
     * @param numRows - the number of rows in the submatrix
     * @param numCols - the number of columns in the submatrix
     * @return A new Matrix object representing the submatrix
     */
    DiagonalMatrix subMatrix(int startRow, int startCol, int numRows, int numCols) const {
        if (startRow < 0 || startRow >= this->getRows() || startCol < 0 || startCol >= this->getCols()) {
            throw std::out_of_range("Invalid start row or column index.");
        }
        if (numRows < 0 || numRows > this->getRows() || numCols < 0 || numCols > this->getCols()) {
            throw std::out_of_range("Invalid number of rows or columns.");
        }
        if (startRow + numRows > this->getRows() || startCol + numCols > this->getCols()) {
            throw std::out_of_range("Submatrix dimensions exceed the original matrix dimensions.");
        }

        DiagonalMatrix result(numRows, T());
        for (int i = 0; i < numRows; ++i) {
            result.setVal(i, i, this->getVal(startRow + i, startCol + i));
        }

        return result;
    }

    /**
     * @brief Calculates the complement of the SquareMatrix
     * @return the complement of the SquareMatrix
     */
    DiagonalMatrix complement() const {
        if (this->det() == T(0)) {
            throw std::invalid_argument("Matrix is singular (determinant is 0)!");
        }

        DiagonalMatrix result(this->getRows(), T(0));

        for (int i = 0; i < this->getRows(); ++i) {
            result.setVal(i, i, this->cofactor(i, i) * ((i % 2 == 0) ? 1 : -1));
        }

        return result;
    }

    /**
     * @brief Calculates the inverse of the SquareMatrix
     * @return the inverse of the SquareMatrix
     */
    DiagonalMatrix inversed() const {
        if (this->det() == T(0)) {
            throw std::invalid_argument("Matrix is singular (determinant is 0)!");
        }

        DiagonalMatrix result(this->getRows());
        for (int i = 0; i < this->getRows(); ++i) {
            result.setVal(i, i, T(1) / this->getVal(i, i));
        }

        return result;
    }

    /**
     * @brief Raises the SquareMatrix to the power of n
     * @param n - the power to raise the SquareMatrix to
     * @return the SquareMatrix raised to the power of n
     */
    DiagonalMatrix power(int n) const {
        if (n < 0) {
            throw std::invalid_argument("Exponent must be non-negative.");
        }

        DiagonalMatrix result(this->getRows(), T(0));
        for (int i = 0; i < this->getRows(); ++i) {
            result.setVal(i, i, std::pow(this->getVal(i, i), n));
        }

        return result;
    }

    
};
}
