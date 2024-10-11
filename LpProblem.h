#ifndef LpProblem_H
#define LpProblem_H

#include <vector>
#include "LpVariable.h"

class SimplexAlgorithm;

class LpProblem {
    friend class SimplexAlgorithm;

    std::vector<double> obj_coefs;
    double obj_b;
    std::vector<std::vector<double>> coefs;
    std::vector<double> b;

    std::vector<LpVariable> variables;
    unsigned int num_variables;

    bool objective_inited;
    void _add_constraint(const Expression& constraint);

    void transform_lower_bound(const int& idx);
public:
    enum ObjectiveType {MAXIMIZE, MINIMIZE};
    ObjectiveType objective;

    LpProblem(const ObjectiveType& obj, const unsigned int& num_vars);
    LpProblem(const ObjectiveType& obj, const Expression& objective_func);
    LpProblem(const LpProblem& other) = default;

    void add_constraint(const Expression& constraint);
    void set_objective(const Expression& objective);
    [[nodiscard]] const std::vector<LpVariable>& get_variables() const;

    void set_variable_bound(const int& id, const double& low, const double& up);
    void set_variable_lower_bound(const int& id, const double& low);
    void set_variable_upper_bound(const int& id, const double& up);

    LpProblem& operator+=(const Expression& constraint);

    void display() const;
};

#endif