//
// Created by QY.
//

#ifndef SimplexAlgorithm_H
#define SimplexAlgorithm_H

#include "LpProblem.h"

class SimplexAlgorithm {

    std::vector<unsigned int> basis;  // The indices of the basic variables
    std::vector<bool> is_basis;  // Keeps track of whether a variable is basic
    std::vector<double> solution;  // Stores the current solution of the LP problem
    LpProblem lp;  // The linear programming problem instance

    bool reversed;  // Flag to indicate if the problem is reversed

    void transform_less_equal(const unsigned int &idx);
    void transform_upper_bound(const int& idx);
    int check_column();
    int check_row(const int& idx);

    void transform_problem();
public:
    explicit SimplexAlgorithm(const LpProblem& problem);

    double solve();
};

#endif