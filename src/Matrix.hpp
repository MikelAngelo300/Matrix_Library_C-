/**
 * @file Matrix.hpp
 * @brief A templated matrix class for general-purpose mathematical operations.
 *
 * This header file defines a generic `Matrix` class that supports arithmetic types (and custom types 
 * that meet certain requirements). The class provides a variety of matrix operations, including:
 * 
 * - Basic arithmetic operations (+, -, *, scalar operations)
 * - Accessors and mutators for elements, rows, and columns
 * - Matrix resizing, transposition, and rotation
 * - Initialization and filling methods (e.g., zeros, ones, random values)
 * - Statistical operations (e.g., sum, mean, min, max)
 * - Normalization and submatrix extraction
 * 
 * The class ensures safety through bounds checking and supports additional functionality like 
 * Gaussian random filling and row/column insertion or deletion.
 * 
 * @tparam T The data type of the matrix elements (must be arithmetic or satisfy custom requirements).
 * 
 * @note Users must manage the header file `ComplexNumber.hpp` dependency for complex number support.
 */


#include "../lib/ComplexNumber.hpp"
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <random>

namespace rm {
template <typename T>
requires std::is_arithmetic_v<T>
class Matrix {
private:
    int m_rows, m_cols;
    T** m_data;

public:

    // Constructor wiht initial value
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

    //Contructor with initializer list
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

    //Copy constructor
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

    // Destructor
    ~Matrix() {
        for (int i = 0; i < m_rows; ++i) {
            delete[] m_data[i];
        }
        delete[] m_data;
    }

    // Getter for number of rows
    int getRows() const { return m_rows; }

    // Getter for number of columns
    int getCols() const { return m_cols; }

    // Getter for data
    T** getData() const { return m_data; }

    // Getter for specific value
    T getVal(int i, int j) const {
        if (i < 0 || i >= m_rows || j < 0 || j >= m_cols) {
            throw std::out_of_range("Matrix indices are out of range.");
        }
        return m_data[i][j];
    }

    // Setter for specific value
    void setVal(int i, int j, T newVal) {
        if (i < 0 || i >= m_rows || j < 0 || j >= m_cols) {
            throw std::out_of_range("Matrix indices are out of range.");
        }
        m_data[i][j] = newVal;
    }

    // Fill entire matrix with a specific value
    virtual void fill(T value) {
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] = value;
            }
        }
    }

    // Fill with random values within a specified range
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

    //Fill the matrix with Gaussian distribution numbers numbers 
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

    //Fill the matrix with 0
    virtual void zeros() {
        this->fill(0);
    }

    //Fill the matrix with 1
    virtual void ones() {
        this->fill(1);
    }

    //Calculate sum of every element in matrix
    T sum() const {
        T sum=0;
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                sum += m_data[i][j];
            }
        }
        return sum;
    }

    //Calculate mean value of every element in matrix
    T mean() const {
        return (this->sum()/(m_cols*m_rows));
    }

    //Find minimal value in matrix
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

    //Find maximal value in matrix
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

    //Normalize matrix (rescale it with numbers between 0 and 1)
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


    //returns size of matrix in pair <rows,cols>
    std::pair<int, int> size() const {
        return {m_rows, m_cols};
    }

    //Extracts submatrix  
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

    //Adds row to the matrix
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

    //Adds column to the matrix
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

    //Removes row from the matrix
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

    //Removes column form the matrix
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

    //Multiplies matrix by a scalar
    void scale(T scalar) {
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] *= scalar;
            }
        }
    }

    // Transposes a matrix
    Matrix transpose() const {
        Matrix<T> transposed(m_cols, m_rows); 

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                transposed.m_data[j][i] = m_data[i][j];
            }
        }

        return transposed;
    }

    // Transposes a matrix to the same one
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
    
    //Rotates matrix by 90 degrees clockwise
    Matrix rotate90() const {
        Matrix result(m_cols, m_rows);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                result.m_data[j][m_rows - 1 - i] = m_data[i][j];
            }
        }

        return result;
    }
    // Rotates matrix by 180 degrees clockwise
    Matrix rotate180() const {
        return this->rotate90().rotate90();
    }

    // Rotates matrix 270 degrees clockwise
    Matrix rotate270() const {
        return this->rotate180().rotate90();
    }

    // Matrix addition
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

    //Matrix addition element wise
    Matrix operator+(const T other) const {
        
        Matrix result(m_rows, m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                result.m_data[i][j] = m_data[i][j] + other;
            }
        }

        return result;
    }

    // Matrix addition assignment
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

    // Matrix addition assignment element wise
    Matrix& operator+=(const T other) {

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] += other;
            }
        }

        return *this;
    
    }

    // Matrix subtraction
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

    // Matrix substraction element wise
    Matrix operator-(const T other) const {
        
        Matrix result(m_rows, m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                result.m_data[i][j] = m_data[i][j] - other;
            }
        }

        return result;
    }

    // Matrix subtraction assignment
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

    // Matrix subtraction assignment element wise
    Matrix& operator-=(const T other) {

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] -= other;
            }
        }

        return *this;
    
    }

    // Matrix multiplication
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

    // Matrix multiplication element wise
    Matrix operator*(const T other) const {
        
        Matrix result(m_rows, other.m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < other.m_cols; ++j) {
                result.m_data[i][j] = m_data[i][j] * other;
            }
        }

        return result;
    }

    // Matrix multiplication assignment
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

    // Matrix multiplication assignment element wise
    Matrix& operator*=(const T other) {
        
        Matrix result(m_rows, other.m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < other.m_cols; ++j) {
                result.m_data[i][j] = m_data[i][j] * other;
            }
        }

        *this = result;

        return *this;
    }

    // Matrix division element wise
    Matrix operator/(const T other) const {
        
        Matrix result(m_rows, other.m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < other.m_cols; ++j) {
                result.m_data[i][j] = m_data[i][j] / other;
            }
        }

        return result;
    }

    // Matrix division assignment elemnt wise
    Matrix& operator/=(const T other) {
        
        Matrix result(m_rows, other.m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < other.m_cols; ++j) {
                result.m_data[i][j] = m_data[i][j] / other;
            }
        }

        *this = result;

        return *this;
    }

    // Assignment operator
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

    // Equality operator
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

    // Inequality operator
    bool operator!=(const Matrix& other) const {
        return !(*this == other);
    }


    // Override of th "<<" operator
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
};
