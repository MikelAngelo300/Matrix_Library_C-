# Matrix Library

This is a C++ library that provides an implementation for matrix operations, including basic matrix manipulations, matrix arithmetic, and operations on square and diagonal matrices. The library supports matrices of numeric types (i.e., arithmetic types).

## Features

- **Matrix Operations**:
  - Matrix addition (`operator+`)
  - Matrix scalar addition (`operator+`)
  - Matrix multiplication
  - Matrix subtraction
  - Matrix transpose
  - Matrix scaling (multiplying all elements by a scalar)
  
- **Square Matrix Operations**:
  - Determinant calculation
  - Cofactor and cofactor matrix
  - Matrix inverse (if determinant is non-zero)
  - Matrix power (exponentiation of matrix)
  - Trace of the matrix (sum of diagonal elements)

- **Diagonal Matrix**:
  - Diagonal matrix creation
  - Access and modification of diagonal elements
  - Identity matrix check
  
- **Matrix Validation**:
  - Dimension checks
  - Square matrix check
  - Diagonal matrix validation

## Classes

### `Matrix<T>`

The base class for all matrix types. It provides fundamental matrix operations and data storage for a 2D array of elements.

#### Constructors
- `Matrix(int rows, int cols, T initValue)` - Constructs a matrix with specified dimensions and initializes all elements to a given value.
- `Matrix(int rows, int cols, std::initializer_list<std::initializer_list<T>> values)` - Constructs a matrix and initializes it using an initializer list.

#### Operators
- `operator+` (Matrix + Matrix) - Adds two matrices element-wise.
- `operator+` (Matrix + Scalar) - Adds a scalar to every element of the matrix.
- `operator*` - Matrix multiplication.
- `operator-` - Matrix subtraction.

#### Methods
- `getVal(int row, int col)` - Access an element of the matrix at the specified position.
- `setVal(int row, int col, T value)` - Set the value of an element in the matrix.
- `transpose()` - Returns a new matrix that is the transpose of the original matrix.
- `scale(T scalar)` - Scales all elements of the matrix by the given scalar.
- `det()` - Returns the determinant of the matrix.
- `trace()` - Returns the trace of the matrix (sum of diagonal elements).

### `SquareMatrix<T>`

This class inherits from `Matrix<T>` and provides additional methods for square matrices.

#### Constructors
- `SquareMatrix(int dim, T initValue)` - Constructs a square matrix with specified dimensions and initializes all elements to a given value.
- `SquareMatrix(int dim, std::initializer_list<std::initializer_list<T>> values)` - Constructs a square matrix and initializes it using an initializer list.

#### Methods
- `cofactor(int row, int col)` - Returns the cofactor of the matrix at the specified position.
- `complement()` - Returns the cofactor matrix.
- `inversed()` - Returns the inverse of the matrix (if the determinant is non-zero).
- `power(int n)` - Returns the matrix raised to the power `n`.

### `DiagonalMatrix<T>`

This class inherits from `SquareMatrix<T>` and is specialized for diagonal matrices.

#### Constructors
- `DiagonalMatrix(int dim, T initValue)` - Creates a diagonal matrix with the specified dimension, with all diagonal elements initialized to a given value.
- `DiagonalMatrix(int dim, std::initializer_list<T> values)` - Creates a diagonal matrix with the specified dimension, with the diagonal elements specified in the initializer list.
- `DiagonalMatrix(const SquareMatrix<T>& other)` - Creates a diagonal matrix from an existing square matrix.

#### Methods
- `getDiagonal()` - Returns the diagonal elements of the matrix.
- `setDiagonal(std::initializer_list<T> values)` - Sets the diagonal elements of the matrix.
- `isIdentity()` - Checks if the matrix is an identity matrix.

## Requirements

- C++20 or higher.
- A C++ compiler such as GCC, Clang, or MSVC.

## Compilation and Usage

To compile the library:

```bash
g++ -std=c++20 -o matrix_program main.cpp
./matrix_program
