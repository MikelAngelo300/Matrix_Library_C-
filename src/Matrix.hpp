#include <iostream>
#include <stdexcept>

template<typename T>
class Matrix {
private:
    int rows, cols;
    T** m_data; // Pointer to the dynamically allocated 2D array

public:
    // Constructor
    Matrix(int rows, int cols, T initValue = T())
        : rows(rows), cols(cols) {
        // Dynamic memory allocation
        m_data = new T*[rows];
        for (int i = 0; i < rows; ++i) {
            m_data[i] = new T[cols];
            for (int j = 0; j < cols; ++j) {
                m_data[i][j] = initValue; // Initialize with default value
            }
        }
    }

    // Destructor (frees memory)
    ~Matrix() {
        for (int i = 0; i < rows; ++i) {
            delete[] m_data[i];
        }
        delete[] m_data;
    }

    int getRows() const { return rows; }

    // Getter for columns
    int getCols() const { return cols; }

    // Getter for data
    T** getData() const { return m_data; }

    // Getter for value
    T getVal(int i, int j) const {
        // Check if indices are within range
        if (i >= rows || j >= cols || i < 0 || j < 0)
            std::cout << "Out of range" << std::endl;
        return m_data[i][j];  // Changed to type T
    }

    // Setter for value
    void setVal(int i, int j, T newVal) {
        if (i >= rows || j >= cols || i < 0 || j < 0)
            std::cout << "Out of range" << std::endl;
        m_data[i][j] = newVal;  // Changed to type T
    }

    template<typename U> // Declaration as a template
    friend std::ostream& operator<<(std::ostream& os, const Matrix<U>& M);

};

// Operator << for printing the matrix
template<typename U>
std::ostream& operator<<(std::ostream& os, const Matrix<U>& M) {
    for (int i = 0; i < M.getRows(); ++i) {      // Loop through rows
        for (int j = 0; j < M.getCols(); ++j) {  // Loop through columns
            os << M.getData()[i][j] << " ";      // Displaying element
        }
        os << std::endl;  // New line after each row
    }
    return os;  // Return the stream to allow further operations on the stream
}
