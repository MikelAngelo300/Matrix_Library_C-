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

        m_data = new T*[rows];   // Dynamic memory allocation

        for (int i = 0; i < rows; ++i) {
            m_data[i] = new T[cols];
            for (int j = 0; j < cols; ++j) {
                m_data[i][j] = initValue; // Initialize with default value
            }
        }
    }

    ~Matrix() {      // Destructor (frees memory)
        for (int i = 0; i < rows; ++i) {
            delete[] m_data[i];
        }
        delete[] m_data;
    }

    int getRows() const { return rows; }     // Getter for rows

    int getCols() const { return cols; }     // Getter for columns

    T** getData() const { return m_data; }   // Getter for data

    T getVal(int i, int j) const {           // Getter for value
        if (i >= rows || j >= cols || i < 0 || j < 0) // Check if indices are within range
            std::cout << "Out of range" << std::endl;
        return m_data[i][j];  // Changed to type T
    }

    void setVal(int i, int j, T newVal) {    // Setter for value
        if (i >= rows || j >= cols || i < 0 || j < 0)
            std::cout << "Out of range" << std::endl;
        m_data[i][j] = newVal;  // Changed to type T
    }

    template<typename U>
    friend std::ostream& operator<<(std::ostream& os, const Matrix<U>& M);

};


template<typename U>
std::ostream& operator<<(std::ostream& os, const Matrix<U>& M) {  // Operator << for printing the matrix
    for (int i = 0; i < M.getRows(); ++i) {      // Loop through rows
        for (int j = 0; j < M.getCols(); ++j) {  // Loop through columns
            os << M.getData()[i][j] << " ";      // Displaying element
        }
        os << std::endl;  // New line after each row
    }
    return os;  // Return the stream to allow further operations on the stream
}
