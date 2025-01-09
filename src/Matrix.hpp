/**
 * @file Matrix.hpp
 * @author Ryszard Mleczko
 * @date Jan 2025
 * @brief A templated matrix class for general-purpose mathematical operations.
 */


#include "../lib/ComplexNumber.hpp"
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <random>

namespace rm {
template <typename T>
requires std::is_arithmetic_v<T>
/**
 * @class Matrix
 * @brief Represents a Matrix with various methods
 */
class Matrix {
private:
    /** Number of rows */
    int m_rows;
    /** Number of columns */
    int m_cols;
    /** Data with which matrix is filled */
    T** m_data;

public:

    /**
     * @brief Parametric constructor of Matrix with default value
     * @param rows - number of rows with which Matrix is initialized
     * @param cols - number of columns with which Matrix is initialized
     * @param initValue - default value with which Matrix is initialized
     */
    Matrix(int rows, int cols, T initValue = T())
        : m_rows(rows), m_cols(cols) {
        if (m_rows <= 0 || m_cols <= 0) {
            throw std::invalid_argument("Matrix dimensions must be positive.");
        }

        m_data = new T*[m_rows];
        for (int i = 0; i < m_rows; ++i) {
            m_data[i] = new T[m_cols];
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] = initValue;
            }
        }
    }

    /**
     * @brief Parametric constructor of Matrix with initializer list
     * @param rows - number of rows with which Matrix is initialized
     * @param cols - number of columns with which Matrix is initialized
     * @param values - values with which Matrix is initialized
     */
    Matrix(int rows, int cols, std::initializer_list<std::initializer_list<T>> values)
        : m_rows(rows), m_cols(cols) {
        if (m_rows <= 0 || m_cols <= 0) {
            throw std::invalid_argument("Matrix dimensions must be positive.");
        }
        if (values.size() != m_rows || values.begin()->size() != m_cols) {
            throw std::invalid_argument("Initializer list dimensions must match the matrix dimensions.");
        }

        m_data = new T*[m_rows];
        auto row_it = values.begin();
        for (int i = 0; i < m_rows; ++i) {
            m_data[i] = new T[m_cols];
            auto col_it = row_it->begin();
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] = *col_it++;
            }
            ++row_it;
        }
    }

    /**
     * @brief Copy constructor of Matrix
     * @param other - Matrix of which is copy being made of 
     */
    Matrix(const Matrix& other)
        : m_rows(other.m_rows), m_cols(other.m_cols) {
        if (m_rows <= 0 || m_cols <= 0) {
            throw std::invalid_argument("Matrix dimensions must be positive.");
        }
        m_data = new T*[m_rows];
        for (int i = 0; i < m_rows; ++i) {
            m_data[i] = new T[m_cols];
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] = other.m_data[i][j];
            }
        }
    }

    /**
     * @brief Default destructor 
     */
    ~Matrix() {
        for (int i = 0; i < m_rows; ++i) {
            delete[] m_data[i];
        }
        delete[] m_data;
    }

    /**
     * Getters
     */

    /**
     * @brief Gets number of rows
     * @return number of rows
     */
    int getRows() const { return m_rows; }

    /**
     * @brief Gets number of columns
     * @return number of cloumns
     */
    int getCols() const { return m_cols; }

    /**
     * @brief Gets all data
     * @return data
     */
    T** getData() const { return m_data; }

    /**
     * @brief Gets data from a specific location
     * @return data from specific location
     */
    T getVal(int i, int j) const {
        if (i < 0 || i >= m_rows || j < 0 || j >= m_cols) {
            throw std::out_of_range("Matrix indices are out of range.");
        }
        return m_data[i][j];
    }

    /**
     * @brief Gets size of matrix in pair <rows,cols>
     * @return size of matrix
     */
    std::pair<int, int> getSize() const {
        return {m_rows, m_cols};
    }

    /**
     * Setters
     */

    /**
     * @brief Sets value at specific location
     * @param i - row number 
     * @param j - column number
     * @param newVal - new Value that is being set 
     */
    void setVal(int i, int j, T newVal) {
        if (i < 0 || i >= m_rows || j < 0 || j >= m_cols) {
            throw std::out_of_range("Matrix indices are out of range.");
        }
        m_data[i][j] = newVal;
    }

    /**
     * Fillers
     */

    /**
     * @brief Fill entire matrix with a specific value
     * @param value - value that is being set into every location 
     */
    virtual void fill(T value) {
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] = value;
            }
        }
    }

    /**
     * @brief Fill entire matrix with zeros 
     */
    virtual void zeros() {
        this->fill(0);
    }

    /**
     * @brief Fill entire matrix with ones 
     */
    virtual void ones() {
        this->fill(1);
    }

    /**
     * @brief Fill with random values within a specified range
     * @param min - minimal range value 
     * @param max - maximal range value
     */
    virtual void randomFill(int min, int max) {
        if (min > max) {
            throw std::invalid_argument("min cannot be greater than max.");
        }

        static bool seeded = false;
        if (!seeded) {
            srand(static_cast<unsigned int>(time(0)));
            seeded = true;
        }

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] = rand() % (max - min + 1) + min;
            }
        }
    }

    /**
     * @brief Fill the matrix with Gaussian distribution numbers numbers
     * @param mean - mean of the Gaussian distribution
     * @param stddev - standard deviation of the Gaussian distribution
     */
    virtual void randomFillGaussian(T mean, T stddev) {
        static_assert(std::is_floating_point<T>::value, "T must be a floating-point type.");
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<T> dist(mean, stddev);


        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] = dist(gen);
            }
        }
    }

    /**
     * Arithmetic methods
     */

    /**
     * @brief Calculates sum of every element in Matrix
     * @return sum
     */
    T sum() const {
        T sum=0;
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                sum += m_data[i][j];
            }
        }
        return sum;
    }

    /**
     * @brief Calculates mean value of every element in Matrix
     * @return mean value
     */
    T mean() const {
        return (this->sum()/(m_cols*m_rows));
    }

    /**
     * @brief Finds element with minimal value in Matrix
     * @return minimal value
     */
    T min() const {
        if (m_rows == 0 || m_cols == 0) {
            throw std::runtime_error("Matrix is empty, cannot determine the minimum value.");
        }

        T buf = m_data[0][0];

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                if (buf>m_data[i][j]){
                    buf=m_data[i][j];
                }
            }
        }
        return buf;
    }

    /**
     * @brief Finds element with maximal value in Matrix
     * @return maximal value
     */
    T max() const {
        if (m_rows == 0 || m_cols == 0) {
            throw std::runtime_error("Matrix is empty, cannot determine the minimum value.");
        }

        T buf = m_data[0][0];

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                if (buf<m_data[i][j]){
                    buf=m_data[i][j];
                }
            }
        }
        return buf;
    }

    /**
     * @brief Normalizes Matrix (rescales it with numbers between 0 and 1)
     */
    void normalize() const { 
        T minValue = this->min();
        T maxValue = this->max();

        if (maxValue == minValue) {
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] /= maxValue;
            }
        } 
        }

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] = (m_data[i][j] - minValue) / (maxValue - minValue);
            }
        }
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
    Matrix subMatrix(int startRow, int startCol, int numRows, int numCols) const {
        if (startRow < 0 || startCol < 0 || startRow + numRows > m_rows || startCol + numCols > m_cols) {
            throw std::out_of_range("Submatrix dimensions exceed matrix boundaries.");
        }
        
        Matrix result(numRows, numCols);
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                result.m_data[i][j] = m_data[startRow + i][startCol + j];
            }
        }

        return result;
    }

    /** 
     * @brief Multiplies Matrix by a scalar
     * @param scalar - number by which every element is being multiplied
     */
    void scale(T scalar) {
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] *= scalar;
            }
        }
    }

    /** 
     * @brief Transposes the current matrix 
     * @return A new Matrix object representing the transposed matrix
     */
    Matrix transpose() const {
        Matrix<T> transposed(m_cols, m_rows); 

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                transposed.m_data[j][i] = m_data[i][j];
            }
        }

        return transposed;
    }

    /**
     * @brief Transposes the matrix in place, modifying the original matrix 
     */
    void selfTranspose() {
        T** newData = new T*[m_cols];
        for (int i = 0; i < m_cols; ++i) {
            newData[i] = new T[m_rows];
            for (int j = 0; j < m_rows; ++j) {
                newData[i][j] = m_data[j][i];
            }
        }

        for (int i = 0; i < m_rows; ++i) {
            delete[] m_data[i];
        }
        delete[] m_data;
        std::swap(m_rows, m_cols);
        m_data = newData;
    }

    /**
     * Size manipulation methods
     */
    
    /**
     * @brief Adds a new row to the matrix at the specified position
     * @param position - the index at which to add the new row
     * @param newRow - the new row to be added, provided as an initializer list 
     */
    void addRow(int position, std::initializer_list<T> newRow) {
        if (position < 0 || position > m_rows) {
            throw std::out_of_range("Invalid row position.");
        }

        if (newRow.size() > m_cols) {
            throw std::invalid_argument("New row has more columns than the matrix.");
        }

        T** newData = new T*[m_rows + 1];
        for (int i = 0; i < m_rows + 1; ++i) {
            newData[i] = new T[m_cols];
        }

        for (int i = 0; i < position; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                newData[i][j] = m_data[i][j];
            }
        }

        for (int i = position; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                newData[i + 1][j] = m_data[i][j];
            }
        }

        int col = 0;
        for (const auto& value : newRow) {
            newData[position][col++] = value;
        }
        for (; col < m_cols; ++col) {
            newData[position][col] = T();
        }

        for (int i = 0; i < m_rows; ++i) {
            delete[] m_data[i];
        }
        delete[] m_data;

        m_data = newData;
        ++m_rows; 
    }

    /**
     * @brief Adds a new column to the matrix at the specified position
     * @param position - the index at which to add the new column
     * @param newCol - the new column to be added, provided as an initializer list 
     */
    void addCol(int position, std::initializer_list<T> newCol) {
        if (position < 0 || position > m_cols) {
            throw std::out_of_range("Invalid column position.");
        }

        if (newCol.size() > m_rows) {
            throw std::invalid_argument("New column has more rows than the matrix.");
        }

        T** newData = new T*[m_rows];
        for (int i = 0; i < m_rows; ++i) {
            newData[i] = new T[m_cols + 1];
        }

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < position; ++j) {
                newData[i][j] = m_data[i][j];
            }
            for (int j = position; j < m_cols; ++j) {
                newData[i][j + 1] = m_data[i][j];
            }
        }

        int row = 0;
        for (const auto& value : newCol) {
            newData[row++][position] = value;
        }
        for (; row < m_rows; ++row) {
            newData[row][position] = T();
        }

        for (int i = 0; i < m_rows; ++i) {
            delete[] m_data[i];
        }
        delete[] m_data;

        m_data = newData;
        ++m_cols;
    }

    /**
     * @brief Removes a row from the matrix at the specified position
     * @param rowIndex - the index at which to remove the row 
     */
    void removeRow(int rowIndex) {
        if (rowIndex < 0 || rowIndex >= m_rows) {
            throw std::out_of_range("Invalid row index.");
        }

        for (int i = rowIndex; i < m_rows - 1; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] = m_data[i + 1][j];
            }
        }

        for (int j = 0; j < m_cols; ++j) {
            m_data[m_rows - 1][j] = T(); 
        }
        
        --m_rows;
    }

    /**
     * @brief Removes a columne from the matrix at the specified position
     * @param colIndex - the index at which to remove the column
     */
    void removeColumn(int colIndex) {
        if (colIndex < 0 || colIndex >= m_cols) {
            throw std::out_of_range("Invalid column index");
        }

        for (int i = 0; i < m_rows; ++i) {
            for (int j = colIndex; j < m_cols - 1; ++j) {
                m_data[i][j] = m_data[i][j + 1];
            }
        }

        for (int i = 0; i < m_rows; ++i) {
            m_data[i][m_cols - 1] = T();
        }

        --m_cols;
    }

    /**
     * Matrix rotations
     */

    /**
     * @brief Rotates matrix by 90 degrees clockwise 
     */
    Matrix rotate90() const {
        Matrix result(m_cols, m_rows);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                result.m_data[j][m_rows - 1 - i] = m_data[i][j];
            }
        }

        return result;
    }

    /**
     * @brief Rotates matrix by 180 degrees clockwise 
     */
    Matrix rotate180() const {
        return this->rotate90().rotate90();
    }

    /**
     * @brief Rotates matrix by 270 degrees clockwise 
     */
    Matrix rotate270() const {
        return this->rotate180().rotate90();
    }

    /**
     * Operators
     */
    
    /** 
     * @brief Adds two matrices element-wise
     * @param other - the matrix to add to the current matrix
     * @return A new Matrix object representing the sum of the two matrices
     */
    Matrix operator+(const Matrix& other) const {
        if (m_rows != other.m_rows || m_cols != other.m_cols) {
            throw std::invalid_argument("Matrix dimensions must be equal.");
        }

        Matrix result(m_rows, m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                result.m_data[i][j] = m_data[i][j] + other.m_data[i][j];
            }
        }

        return result;
    }

    /**
     * @brief Adds a scalar value to each element of the matrix
     * @param other - the scalar value to be added to each element
     * @return A new Matrix object with the scalar added to each element
     */
    Matrix operator+(const T other) const {
        
        Matrix result(m_rows, m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                result.m_data[i][j] = m_data[i][j] + other;
            }
        }

        return result;
    }

    /** 
     * @brief Adds another matrix to the current matrix element-wise
     * @param other - the matrix to add to the current matrix
     * @return A reference to the current matrix after addition
     */
    Matrix& operator+=(const Matrix& other) {
        if (m_rows != other.m_rows || m_cols != other.m_cols) {
            throw std::invalid_argument("Matrix dimensions must be equal.");
        }

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] += other.m_data[i][j];
            }
        }

        return *this;
    }

  /**
    * @brief Adds a scalar value to each element of the current matrix
    * @param other - the scalar value to be added to each element
    * @return A reference to the current matrix after addition
    */
    Matrix& operator+=(const T other) {

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] += other;
            }
        }

        return *this;
    
    }

    /** 
     * @brief Subtracts two matrices element-wise
     * @param other - the matrix to subtract from the current matrix
     * @return A new Matrix object representing the subtraction of the two matrices
     */
    Matrix operator-(const Matrix& other) const {
        if (m_rows != other.m_rows || m_cols != other.m_cols) {
            throw std::invalid_argument("Matrix dimensions must be equal.");
        }

        Matrix result(m_rows, m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                result.m_data[i][j] = m_data[i][j] - other.m_data[i][j];
            }
        }

        return result;
    }

    /**
     * @brief Subtracts a scalar value from each element of the matrix
     * @param other - the scalar value to be subtracted from each element
     * @return A new Matrix object with the scalar subtracted from each element
     */
    Matrix operator-(const T other) const {
        
        Matrix result(m_rows, m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                result.m_data[i][j] = m_data[i][j] - other;
            }
        }

        return result;
    }

    /**
     * @brief Adds another matrix to the current matrix element-wise
     * @param other - the matrix to add to the current matrix
     * @return A reference to the current matrix after addition   
     */
    Matrix& operator-=(const Matrix& other) {
        if (m_rows != other.m_rows || m_cols != other.m_cols) {
            throw std::invalid_argument("Matrix dimensions must be equal.");
        }

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] -= other.m_data[i][j];
            }
        }

        return *this;
    }

    /**
     * @brief Adds a scalar value to each element of the current matrix
     * @param other - the scalar value to be added to each element
     * @return A reference to the current matrix after addition
     */
    Matrix& operator-=(const T other) {

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] -= other;
            }
        }

        return *this;
    }

    /**
     * @brief Multiplies the current matrix with another matrix
     * @param other - the matrix to multiply with
     * @return A new Matrix object representing the product of the two matrices
     * @throws std::invalid_argument if the dimensions do not allow multiplication
     */
    Matrix operator*(const Matrix& other) const {
        if (m_cols != other.m_rows) {
            throw std::invalid_argument("Matrix dimensions do not allow multiplication.");
        }

        Matrix result(m_rows, other.m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < other.m_cols; ++j) {
                T sum = 0;
                for (int k = 0; k < m_cols; ++k) {
                    sum += m_data[i][k] * other.m_data[k][j];
                }
                result.m_data[i][j] = sum;
            }
        }

        return result;
    }

    /**
     * @brief Multiplies each element of the matrix by a scalar value
     * @param other - the scalar value to multiply by
     * @return A new Matrix object with each element multiplied by the scalar
     */
    Matrix operator*(const T other) const {
        Matrix result(m_rows, m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                result.m_data[i][j] = m_data[i][j] * other;
            }
        }

        return result;
    }

    /**
     * @brief Multiplies the current matrix with another matrix and assigns the result to the current matrix
     * @param other - the matrix to multiply with
     * @return A reference to the current matrix after multiplication
     * @throws std::invalid_argument if the dimensions do not allow multiplication
     */
    Matrix& operator*=(const Matrix& other) {
        if (m_cols != other.m_rows) {
            throw std::invalid_argument("Matrix dimensions do not allow multiplication.");
        }

        Matrix result(m_rows, other.m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < other.m_cols; ++j) {
                T sum = 0;
                for (int k = 0; k < m_cols; ++k) {
                    sum += m_data[i][k] * other.m_data[k][j];
                }
                result.m_data[i][j] = sum;
            }
        }

        *this = result;

        return *this;
    }

    /**
     * @brief Multiplies each element of the current matrix by a scalar value and assigns the result to the current matrix
     * @param other - the scalar value to multiply by
     * @return A reference to the current matrix after multiplication
     */
    Matrix& operator*=(const T other) {
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] *= other;
            }
        }

        return *this;
    }

    /**
     * @brief Divides each element of the matrix by a scalar value
     * @param other - the scalar value to divide by
     * @return A new Matrix object with each element divided by the scalar
     */
    Matrix operator/(const T other) const {
        Matrix result(m_rows, m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                result.m_data[i][j] = m_data[i][j] / other;
            }
        }

        return result;
    }

    /**
     * @brief Divides each element of the current matrix by a scalar value and assigns the result to the current matrix
     * @param other - the scalar value to divide by
     * @return A reference to the current matrix after division
     */
    Matrix& operator/=(const T other) {
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] /= other;
            }
        }

        return *this;
    }

    /**
     * @brief Assigns the values from another matrix to the current matrix
     * @param other - the matrix to assign from
     * @return A reference to the current matrix after assignment
     */
    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            // Deallocate existing memory
            for (int i = 0; i < m_rows; ++i) {
                delete[] m_data[i];
            }
            delete[] m_data;

            m_rows = other.m_rows;
            m_cols = other.m_cols;
            m_data = new T*[m_rows];
            for (int i = 0; i < m_rows; ++i) {
                m_data[i] = new T[m_cols];
                for (int j = 0; j < m_cols; ++j) {
                    m_data[i][j] = other.m_data[i][j];
                }
            }
        }
        return *this;
    }

    /**
     * @brief Checks if two matrices are equal
     * @param other - the matrix to compare with
     * @return True if the matrices are equal, false otherwise
     */
    bool operator==(const Matrix& other) const {
        if (m_rows != other.m_rows || m_cols != other.m_cols) {
            return false;
        }

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                if (m_data[i][j] != other.m_data[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @brief Checks if two matrices are not equal
     * @param other - the matrix to compare with
     * @return True if the matrices are not equal, false otherwise
     */
    bool operator!=(const Matrix& other) const {
        return !(*this == other);
    }

    /**
     * @brief Outputs the matrix elements to a stream
     * @param os - the output stream
     * @param M - the matrix to output
     * @return The output stream with the matrix elements
     */
    template<typename U>
    friend std::ostream& operator<<(std::ostream& os, const Matrix<U>& M);
};
// Stream output implementation
template<typename U>
std::ostream& operator<<(std::ostream& os, const Matrix<U>& M) {
    for (int i = 0; i < M.getRows(); ++i) {
        for (int j = 0; j < M.getCols(); ++j) {
            os << M.getData()[i][j] << " ";
        }
        os << std::endl;
    }
    return os;
}
}

