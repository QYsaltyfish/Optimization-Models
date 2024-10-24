//
// Created by QY.
//

#ifndef LpProblem_H
#define LpProblem_H

#include <vector>
#include "LpVariable.h"

class SimplexAlgorithm;

class LpProblem {
    friend class SimplexAlgorithm;

    std::vector<double> obj_coefs;  // Coefficients of the objective function
    double obj_b;  // Right-hand side value of the objective function
    std::vector<std::vector<double>> coefs;  // Coefficients of the constraints
    std::vector<double> b;  // Right-hand side values of the constraints

    std::vector<LpVariable> variables;  // List of LP variables
    unsigned int num_variables;  // Number of variables in the LP problem

    bool objective_inited;  // Flag to check if the objective function is initialized
    void _add_constraint(const Expression& constraint);

    void transform_lower_bound(const int& idx);
public:
    enum ObjectiveType {MAXIMIZE, MINIMIZE};
    ObjectiveType objective;

    LpProblem(const ObjectiveType& obj, const unsigned int& num_vars);
    LpProblem(const ObjectiveType& obj, const Expression& objective_func);
    explicit LpProblem(const std::string& csvFilePath);
    LpProblem(const LpProblem& other) = default;

    void add_constraint(const Expression& constraint);
    void set_objective(const Expression& objective);
    [[nodiscard]] const std::vector<LpVariable>& get_variables() const;

    void set_variable_bound(const int& id, const double& low, const double& up);
    void set_bounds(const double& low, const double& up);
    void set_variable_lower_bound(const int& id, const double& low);
    void set_lower_bounds(const double& low);
    void set_variable_upper_bound(const int& id, const double& up);
    void set_upper_bounds(const double& up);

    LpProblem& operator+=(const Expression& constraint);

    void display() const;
};

#endif