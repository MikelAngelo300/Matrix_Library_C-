#include "../lib/ComplexNumber.hpp"
#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <ctime>

namespace rm {
    template<typename T>
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
