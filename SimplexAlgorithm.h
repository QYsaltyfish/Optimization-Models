#ifndef SimplexAlgorithm_H
#define SimplexAlgorithm_H

#include "LpProblem.h"

class SimplexAlgorithm {

    std::vector<unsigned int> basis;
    std::vector<bool> is_basis;
    std::vector<double> solution;
    LpProblem lp;

    bool reversed;

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