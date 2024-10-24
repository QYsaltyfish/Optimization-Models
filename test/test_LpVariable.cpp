//
// Created by QY.
//

#include "../LpVariable.h"

int main() {
    LpVariable vars[10];

    std::cout << std::endl << "Generate 10 LP variables:" << std::endl;
    for (const auto& var : vars)
        std::cout << var << ' ';
    std::cout << std::endl << std::endl;

    std::cout << "A simple linear combination of LP variables:" << std::endl;
    Expression expr1 = vars[0] * 2 - vars[1] * 1.2 + 2;
    const Expression expr2 = (vars[1] * 1.2 - 2) * 0.5;
    std::cout << "Expr1: " << expr1 << std::endl;
    std::cout << "Expr2: " << expr2 << std::endl << std::endl;

    std::cout << "Expr1 + x2: " << expr1 + vars[2] << std::endl;
    std::cout << "Expr1 - x2: " << expr1 - vars[2] << std::endl;
    std::cout << "Expr1 + Expr2: " << expr1 + expr2 << std::endl;
    std::cout << "Expr1 - Expr2: " << expr1 - expr2 << std::endl;
    std::cout << "Expr1 + 1.1: " << expr1 + 1.1 << std::endl;
    std::cout << "Expr1 - 1.1: " << expr1 - 1.1 << std::endl << std::endl;

    expr1 += vars[2];
    std::cout << "Expr1 += x2: " << expr1 << std::endl;
    expr1 -= vars[2];
    std::cout << "Expr1 -= x2: " << expr1 << std::endl;
    expr1 += expr2;
    std::cout << "Expr1 += Expr2: " << expr1 << std::endl;
    expr1 -= expr2;
    std::cout << "Expr1 -= Expr2: " << expr1 << std::endl;
    expr1 += 1.1;
    std::cout << "Expr1 += 1.1: " << expr1 << std::endl;
    expr1 -= 1.1;
    std::cout << "Expr1 -= 1.1: " << expr1 << std::endl << std::endl;

    std::cout << "Expr1 * 2: " << expr1 * 2 << std::endl;
    std::cout << "Expr1 / 2: " << expr1 / 2 << std::endl;
    expr1 *= 2;
    std::cout << "Expr1 *= 2: " << expr1 << std::endl;
    expr1 /= 2;
    std::cout << "Expr1 /= 2: " << expr1 << std::endl << std::endl;

    std::cout << "Expr1: " << expr1 << std::endl;
    std::cout << "Expr2: " << expr2 << std::endl << std::endl;

    std::cout << "Expr1 <= Expr2: " << (expr1 <= expr2) << std::endl;
    std::cout << "Expr1 <= x2: " << (expr1 <= vars[2]) << std::endl;
    std::cout << "Expr1 <= 1.1: " << (expr1 <= 1.1) << std::endl << std::endl;

    std::cout << "Expr1 >= Expr2: " << (expr1 >= expr2) << std::endl;
    std::cout << "Expr1 >= x2: " << (expr1 >= vars[2]) << std::endl;
    std::cout << "Expr1 >= 1.1: " << (expr1 >= 1.1) << std::endl << std::endl;

    std::cout << "Expr1 == Expr2: " << (expr1 == expr2) << std::endl;
    std::cout << "Expr1 == x2: " << (expr1 == vars[2]) << std::endl;
    std::cout << "Expr1 == 1.1: " << (expr1 == 1.1) << std::endl << std::endl;
}
