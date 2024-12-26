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

    // Constructor
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
    void fill(T value) {
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] = value;
            }
        }
    }

    // Fill with random values within a specified range
    void randomFill(int min, int max) {
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
    void randomFillGaussian(T mean, T stddev) {
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
    void zeros() {
        this->fill(0);
    }

    //Fill the matrix with 1
    void ones() {
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

    //Find minimal value in matrix
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
    void normalize() { 
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


    void resize(int newRows, int newCols, T defaultValue = T()) {
        m_rows += newRows;
        m_cols += newCols;

        for (int i = newRows; i < m_rows; ++i) {
            for (int j = newCols; j < m_cols; ++j) {
                m_data[i][j] = defaultValue;
            }
        }
    }

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

    void addRow(int position, std::initializer_list<T> newRow) {
        if (position < 0 || position > m_rows) {
            throw std::out_of_range("Invalid row position.");
        }

        // Ensure the new row has the correct number of columns
        if (newRow.size() > m_cols) {
            throw std::invalid_argument("New row has more columns than the matrix.");
        }

        // Create a temporary row array of size m_cols
        T fullRow[m_cols];

        // Copy elements from newRow into fullRow
        size_t j = 0;
        for (const auto& value : newRow) {
            fullRow[j++] = value;
        }

        // Fill the rest of the row with default values if necessary
        for (; j < m_cols; ++j) {
            fullRow[j] = T();  // Default constructor of T
        }

        // Insert the row into the matrix at the specified position
        // Shift rows downwards
        for (int i = m_rows; i > position; --i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] = m_data[i - 1][j];
            }
        }

        // Insert the new row at the specified position
        for (int j = 0; j < m_cols; ++j) {
            m_data[position][j] = fullRow[j];
        }

        // Update the row count
        ++m_rows;
    }



    void addCol(int position, std::initializer_list<T> newCol) {
        if (position < 0 || position > m_cols) {
            throw std::out_of_range("Invalid column position.");
        }

        if (newCol.size() > m_rows) {
            throw std::invalid_argument("New column has more rows than the matrix.");
        }

        T fullCol[m_rows];

        size_t j = 0;
        for (const auto& value : newCol) {
            fullCol[j++] = value;
        }

        for (; j < m_rows; ++j) {
            fullCol[j] = T();
        }

        for (int i = m_cols; i > position; --i) {
            for (int j = 0; j < m_rows; ++j) {
                m_data[j][i] = m_data[j][i - 1];
            }
        }

        for (int j = 0; j < m_rows; ++j) {
            m_data[j][position] = fullCol[j];
        }

        ++m_cols;
    }

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


    void removeColumn(int colIndex) {
        if (colIndex < 0 || colIndex >= m_cols) {
            throw std::out_of_range("Invalid column index");
        }

        // Shift columns to the left to fill the gap created by the removed column
        for (int i = 0; i < m_rows; ++i) {
            for (int j = colIndex; j < m_cols - 1; ++j) {
                m_data[i][j] = m_data[i][j + 1];
            }
        }

        // Clear the last column of each row (optional, depending on your implementation)
        for (int i = 0; i < m_rows; ++i) {
            m_data[i][m_cols - 1] = T();  // Default constructor for type T
        }

        // Update the column count
        --m_cols;
    }

    
    void scale(T scalar) {
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                m_data[i][j] *= scalar;
            }
        }
    }

    Matrix transpose() const {
        // Create a new matrix with swapped dimensions
        Matrix<T> transposed(m_cols, m_rows); 

        // Loop over the original matrix and assign values to the transposed matrix
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                transposed.m_data[j][i] = m_data[i][j];  // Swap rows and columns
            }
        }

        return transposed;
    }

    void selfTranspose() {
        // Utworzenie kopii danych, ponieważ transponowanie zmienia wymiary
        T** newData = new T*[m_cols];
        for (int i = 0; i < m_cols; ++i) {
            newData[i] = new T[m_rows];
            for (int j = 0; j < m_rows; ++j) {
                newData[i][j] = m_data[j][i];
            }
        }

        // Dealokacja starej macierzy
        for (int i = 0; i < m_rows; ++i) {
            delete[] m_data[i];
        }
        delete[] m_data;

        // Zamiana wymiarów i przypisanie nowej macierzy
        std::swap(m_rows, m_cols);
        m_data = newData;
    }

    Matrix rotate90() const {
        // Tworzymy wynikową macierz o wymiarach m_cols x m_rows
        Matrix result(m_cols, m_rows);

        // Iteracja przez oryginalną macierz i przekształcanie elementów
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                // Obrót o 90 stopni zgodnie z zasadami
                result.m_data[j][m_rows - 1 - i] = m_data[i][j];
            }
        }

        return result;
    }

    Matrix rotate180() const {
        return this->rotate90().rotate90();
    }

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

     Matrix operator/(const T other) const {
        
        Matrix result(m_rows, other.m_cols);

        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < other.m_cols; ++j) {
                result.m_data[i][j] = m_data[i][j] / other;
            }
        }

        return result;
    }

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

            // Allocate new memory and copy values
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
        return !(*this == other);  // Delegate to the equality operator
    }

    // Stream output operator
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
