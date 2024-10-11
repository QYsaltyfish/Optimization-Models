#ifndef Matrix_H
#define Matrix_H

#include <iostream>
#include <stdexcept>
#include <string>

template <typename T>
class Matrix {
    void throw_dimension_error(const Matrix &other) const;

    static Matrix LU_inv(Matrix& mat);
    static void LU_decomposition(Matrix& L, Matrix& U);
public:
    int row, col;
    T **data;

    Matrix(const int& r, const int& c);
    Matrix(const Matrix& other);
    [[maybe_unused]] Matrix(T **array, int r, int c);
    ~Matrix();

    static Matrix<T> eye(const int& size);
    static Matrix<T> zeros(const int& r, const int& c);

    Matrix operator+(const Matrix& other) const;
    Matrix operator+(const T& scalar) const;
    Matrix& operator+=(const Matrix& other);
    Matrix& operator+=(const T& scalar);

    Matrix operator-(const Matrix& other) const;
    Matrix operator-(const T& scalar) const;
    Matrix& operator-=(const Matrix& other);
    Matrix& operator-=(const T& scalar);

    Matrix operator*(const Matrix& other) const;
    Matrix operator*(const T& scalar) const;
    Matrix& operator*=(const Matrix& other);
    Matrix& operator*=(const T& scalar);

    Matrix operator/(const T& scalar) const;
    Matrix& operator/=(const T& scalar);

    Matrix inv() const;
    Matrix t() const;

    void display() const;
};

template<typename T>
Matrix<T> Matrix<T>::t() const {
    Matrix<T> result = Matrix(col, row);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j)
            result.data[j][i] = data[i][j];
    }

    return result;
}

template<typename T>
Matrix<T>::Matrix(const Matrix &other): row(other.row), col(other.col) {
    data = new T*[row];

    for (int i = 0; i < row; ++i) {
        data[i] = new T[col];
        for (int j = 0; j < col; ++j)
            data[i][j] = other.data[i][j];
    }
}

template<typename T>
void Matrix<T>::LU_decomposition(Matrix &L, Matrix &U) {
    for (int j = 0; j < L.row; ++j) {
        for (int i = j + 1; i < L.row; ++i) {
            T temp = U.data[i][j] / U.data[j][j];
            L.data[i][j] = temp;
            for (int k = j; k < L.row; ++k)
                U.data[i][k] -= temp * U.data[j][k];
        }
    }
}

template<typename T>
Matrix<T> Matrix<T>::LU_inv(Matrix &mat) {
    Matrix<T> L = Matrix::eye(mat.row);
    Matrix<T> result = Matrix::zeros(mat.row, mat.row);

    // Perform LU decomposition to acquire Matrix L and U
    LU_decomposition(L, mat);

    // Solve L * y_i = e_i and U * x_i = y_i
    auto y = new T [mat.row];
    T rhs;
    for (int i = 0; i < mat.row; ++i) {

        // Step 1: Solve L * y_i = e_i
        for (int j = 0; j < mat.row; ++j) {
            rhs = (j == i) ? 1 : 0;
            for (int k = 0; k < j; ++k)
                rhs -= L.data[j][k] * y[k];
            y[j] = rhs;
        }

        // Step 2: Solve U * x_i = y_i
        for (int j = mat.row - 1; j >= 0; --j) {
            for (int k = mat.row - 1; k > j; --k) {
                y[j] -= mat.data[j][k] * y[k];
            }
            y[j] /= mat.data[j][j];
        }

        for (int j = 0; j < mat.row; ++j)
            result.data[j][i] = y[j];
    }

    delete [] y;
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::inv() const {
    if (row != col)
        throw std::invalid_argument("Matrix must be square to compute its inverse.");

    Matrix mat = *this;
    return LU_inv(mat);  // Obtain the inverse matrix using LU decomposition
    // More methods may be added in the future
}

template<typename T>
Matrix<T> Matrix<T>::eye(const int &size) {
    auto result = Matrix<T>(size, size);

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < i; ++j)
            result.data[i][j] = 0;
        for (int j = i + 1; j < size; ++j)
            result.data[i][j] = 0;
        result.data[i][i] = 1;
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::zeros(const int &r, const int &c) {
    auto result = Matrix<T>(r, c);

    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j)
            result.data[i][j] = 0;

    return result;
}

template <typename T>
void Matrix<T>::throw_dimension_error(const Matrix &other) const {
    throw std::invalid_argument(
            "Matrix dimensions do not match, received Matrix 1: (" +
            std::to_string(row) + ", " + std::to_string(col) +
            ") and Matrix 2: (" +
            std::to_string(other.row) + ", " + std::to_string(other.col) +
            ")"
    );
}

template <typename T>
Matrix<T>::Matrix(const int& r, const int& c) {
    row = r;
    col = c;

    data = new T*[row];

    for (int i = 0; i < r; ++i)
        data[i] = new T[col];
}

template <typename T>
[[maybe_unused]] Matrix<T>::Matrix(T** array, int r, int c) {
    row = r;
    col = c;

    data = new T* [r];
    for (int i = 0; i < r; ++i) {
        data[i] = new T [c];
        for (int j = 0; j < c; ++j)
            data[i][j] = array[i][j];
    }
}

template <typename T>
Matrix<T>::~Matrix() {
    for (int i = 0; i < row; ++i) {
        delete[] data[i];
    }
    delete[] data;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix &other) const {
    if (row != other.row || col != other.col)
        throw_dimension_error(other);

    Matrix result(row, col);

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j)
            result.data[i][j] = data[i][j] + other.data[i][j];
    }

    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const T& scalar) const {
    Matrix result(row, col);

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j)
            result.data[i][j] = data[i][j] + scalar;
    }

    return result;
}

template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix &other) {
    if (row != other.row || col != other.col)
        throw_dimension_error(other);

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j)
            data[i][j] += other.data[i][j];
    }

    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator+=(const T& scalar) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j)
            data[i][j] += scalar;
    }

    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix &other) const {
    if (row != other.row || col != other.col)
        throw_dimension_error(other);

    Matrix result(row, col);

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j)
            result.data[i][j] = data[i][j] - other.data[i][j];
    }

    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const T& scalar) const {
    Matrix result(row, col);

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j)
            result.data[i][j] = data[i][j] - scalar;
    }

    return result;
}

template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix &other) {
    if (row != other.row || col != other.col)
        throw_dimension_error(other);

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j)
            data[i][j] -= other.data[i][j];
    }

    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator-=(const T& scalar) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j)
            data[i][j] -= scalar;
    }

    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix &other) const {
    if (col != other.row)
        throw_dimension_error(other);

    Matrix result = Matrix::zeros(row, other.col);

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < other.col; ++j)
            for (int k = 0; k < col; ++k)
                result.data[i][j] += data[i][k] * other.data[k][j];
    }

    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const T& scalar) const {
    Matrix result(row, col);

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j)
            result.data[i][j] = data[i][j] * scalar;
    }

    return result;
}

template <typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix &other) {
    return *this * other;
}

template <typename T>
Matrix<T>& Matrix<T>::operator*=(const T& scalar) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j)
            data[i][j] *= scalar;
    }

    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator/(const T& scalar) const {
    Matrix result(row, col);

    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j)
            result.data[i][j] = data[i][j] / scalar;
    }

    return result;
}

template <typename T>
Matrix<T>& Matrix<T>::operator/=(const T& scalar) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j)
            data[i][j] /= scalar;
    }

    return *this;
}

template <typename T>
void Matrix<T>::display() const {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j)
            std::cout << data[i][j] << ' ';
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

#endif
