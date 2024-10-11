#include "Matrix.h"

int main() {
    Matrix<int> mat1(3, 3);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            mat1.data[i][j] = i * 3 + j;
    }

    std::cout << "Matrix 1:" << std::endl;
    mat1.display();

    const Matrix<int> mat2(3, 3);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            mat2.data[i][j] = i * 3 + j + 1;
    }

    std::cout << "Matrix 2:" << std::endl;
    mat2.display();

    std::cout << "Matrix 1 + Matrix 2:" << std::endl;
    (mat1 + mat2).display();

    std::cout << "Matrix 1 - Matrix 2:" << std::endl;
    (mat1 - mat2).display();

    std::cout << "Matrix 1 + 1:" << std::endl;
    (mat1 + 1).display();

    std::cout << "Matrix 1 - 1:" << std::endl;
    (mat1 - 1).display();

    std::cout << "Matrix 1 += 1:" << std::endl;
    mat1 += 1;
    mat1.display();

    std::cout << "Matrix 1 -= 1:" << std::endl;
    mat1 -= 1;
    mat1.display();

    std::cout << "Matrix 1 += Matrix 2:" << std::endl;
    mat1 += mat2;
    mat1.display();

    std::cout << "Matrix 1 -= Matrix 2:" << std::endl;
    mat1 -= mat2;
    mat1.display();

    std::cout << "Matrix 1 * Matrix 2:" << std::endl;
    (mat1 * mat2).display();

    std::cout << "Matrix 1 * 2:" << std::endl;
    (mat1 * 2).display();

    std::cout << "Matrix 1 *= 2:" << std::endl;
    mat1 *= 2;
    mat1.display();

    std::cout << "Matrix 1 /= 2:" << std::endl;
    mat1 /= 2;
    mat1.display();

    const Matrix<double> mat3(4, 4);
    mat3.data[0][0] = 2; mat3.data[0][1] = 1; mat3.data[0][2] = 1; mat3.data[0][3] = 0;
    mat3.data[1][0] = 4; mat3.data[1][1] = 3; mat3.data[1][2] = 3; mat3.data[1][3] = 1;
    mat3.data[2][0] = 8; mat3.data[2][1] = 7; mat3.data[2][2] = 9; mat3.data[2][3] = 5;
    mat3.data[3][0] = 6; mat3.data[3][1] = 7; mat3.data[3][2] = 9; mat3.data[3][3] = 8;
    std::cout << "Matrix 3:" << std::endl;
    mat3.display();

    std::cout << "The inverse of Matrix 3:" << std::endl;
    (mat3.inv()).display();

    std::cout << "The transpose of Matrix 3:" << std::endl;
    (mat3.t()).display();

    return 0;
}
