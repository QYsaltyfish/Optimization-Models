//
// Created by QY.
//

#include "../TSP.h"
#include <iostream>

int main() {

    // Test Prim algorithm
    std::vector<std::vector<int>> adj_matrix = {
            {0, 23, INT_MAX, INT_MAX, INT_MAX, 28, 36},
            {23, 0, 20, INT_MAX, INT_MAX, INT_MAX, 1},
            {INT_MAX, 20, 0, 15, INT_MAX, INT_MAX, 4},
            {INT_MAX, INT_MAX, 15, 0, 3, INT_MAX, 9},
            {INT_MAX, INT_MAX, INT_MAX, 3, 0, 17, 16},
            {28, INT_MAX, INT_MAX, INT_MAX, 17, 0, 25},
            {36, 1, 4, 9, 16, 25, 0}
    };
    TspProblem problem1 = TspProblem(adj_matrix, 0);
    LKHSolver solver1 = LKHSolver(problem1);
    std::cout << solver1.minimum_spanning_tree();
}