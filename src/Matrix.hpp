#include <iostream>
#include <stdexcept>

template<typename T>
class Matrix {
private:
    int rows, cols;
    T** m_data; // Wskaźnik do dynamicznej tablicy dwuwymiarowej

public:
    // Konstruktor
    Matrix(int rows, int cols, T initValue = T())
        : rows(rows), cols(cols) {
        // Dynamiczna alokacja pamięci
        m_data = new T*[rows];
        for (int i = 0; i < rows; ++i) {
            m_data[i] = new T[cols];
            for (int j = 0; j < cols; ++j) {
                m_data[i][j] = initValue; // Inicjalizujemy wartością domyślną
            }
        }
    }

    // Destruktor (zwalnia pamięć)
    ~Matrix() {
        for (int i = 0; i < rows; ++i) {
            delete[] m_data[i];
        }
        delete[] m_data;
    }

    // Dostęp do elementu
    T& at(int row, int col) {
        if (row >= rows || col >= cols || row < 0 || col < 0)
            std::cout<<"Out of range"<< std::endl;
        return m_data[row][col];
    }

    int getRows() const { return rows; }

    // Getter dla cols
    int getCols() const { return cols; }

    // Getter dla danych
    T** getData() const { return m_data; }

    // Getter dla wartości (stały dostęp)
    T getVal(int i, int j) const {
        // Sprawdzamy, czy indeksy są w zakresie
        if (i >= rows || j >= cols || i < 0 || j < 0)
            std::cout << "Out of range" << std::endl;
        return m_data[i][j];  // Zmieniono na typ T
    }

    // Setter dla wartości
    void setVal(int i, int j, T newVal) {
        if (i >= rows || j >= cols || i < 0 || j < 0)
            std::cout << "Out of range" << std::endl;
        m_data[i][j] = newVal;  // Zmieniono na typ
    }

    template<typename U> // Deklaracja jako szablon
    friend std::ostream& operator<<(std::ostream& os, const Matrix<U>& M);

};

template<typename U>
std::ostream& operator<<(std::ostream& os, const Matrix<U>& M) {
    for (int i = 0; i < M.getRows(); ++i) {      // Pętla po wierszach
        for (int j = 0; j < M.getCols(); ++j) {  // Pętla po kolumnach
            os << M.getData()[i][j] << " ";      // Wyświetlanie elementu
        }
        os << std::endl;  // Nowa linia po każdym wierszu
    }
    return os;  // Zwracamy strumień, aby umożliwić dalsze operacje na strumieniu
}