//
// Created by QY.
//

#include "../TSP.h"
#include <iostream>

int main() {

    // Test Prim algorithm
    std::vector<std::vector<int>> adj_matrix = {
            {0, 23, 10000, 10000, 10000, 28, 36},
            {23, 0, 20, 10000, 10000, 10000, 1},
            {10000, 20, 0, 15, 10000, 10000, 4},
            {10000, 10000, 15, 0, 3, 10000, 9},
            {10000, 10000, 10000, 3, 0, 17, 16},
            {28, 10000, 10000, 10000, 17, 0, 25},
            {36, 1, 4, 9, 16, 25, 0},
    };
    TspProblem problem1 = TspProblem(adj_matrix, 0);
    LKHSolver solver1 = LKHSolver(problem1);
    // std::cout << solver1.ascent();
    solver1.create_candidate_set();

    return 0;
}