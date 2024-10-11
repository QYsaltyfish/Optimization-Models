#include <iostream>
#include "LpProblem.h"

void LpProblem::_add_constraint(const Expression &constraint) {
    int count = 0;
    int id = -1;

    for (int i = 0; i < num_variables; ++i) {
        if (constraint.arr[i] == 0)
            continue;
        ++count;
        id = i;
        if (count == 2)
            break;
    }

    if (count == 2) {
        coefs.push_back(constraint.arr);
        b.push_back(constraint.b);
        if (constraint.exp_type == Expression::EQUAL) {
            auto new_vector = std::vector<double>(num_variables);
            for (int i = 0; i < num_variables; ++i)
                new_vector[i] = -constraint.arr[i];
            coefs.push_back(new_vector);
            b.push_back(-constraint.b);
        }
    }
    else {
        if (id == -1)
            return;
        if (constraint.arr[id] > 0 && constraint.arr[id] < variables[id].upper_bound)
            set_variable_upper_bound(id, constraint.b / constraint.arr[id]);
        else if (constraint.arr[id] < 0 && constraint.arr[id] > variables[id].lower_bound)
            set_variable_lower_bound(id, constraint.b / constraint.arr[id]);
    }
}

void LpProblem::transform_lower_bound(const int &idx) {

    // Case 1: The variable has a lower bound
    if (variables[idx].lower_bound != -std::numeric_limits<double>::infinity()) {

        // Case 1.1: The lower bound is exactly 0
        if (variables[idx].lower_bound == 0)
            return;

        // Case 1.2: The lower bound is not 0
        for (int row = 0; row < coefs.size(); ++row) {
            if (coefs[row][idx] == 0)
                continue;
            b[row] -= coefs[row][idx] * variables[idx].lower_bound;
        }
        obj_b += obj_coefs[idx] * variables[idx].lower_bound;
        variables[idx].lower_bound = 0;
        variables[idx].upper_bound -= variables[idx].lower_bound;
        return;
    }

    // Case 2: The variable doesn't have a lower bound
    // Case 2.1: The variable has an upper bound
    if (variables[idx].upper_bound != std::numeric_limits<double>::infinity()) {
        for (int row = 0; row < coefs.size(); ++row) {
            if (coefs[row][idx] == 0)
                continue;
            b[row] -= coefs[row][idx] * variables[idx].upper_bound;
            coefs[row][idx] = -coefs[row][idx];
        }
        obj_b += obj_coefs[idx] * variables[idx].upper_bound;
        variables[idx].lower_bound = 0;
        variables[idx].upper_bound = std::numeric_limits<double>::infinity();
        return;
    }

    // Case 2.2: The variable doesn't have an upper bound
    ++num_variables;

    obj_coefs.push_back(0);
    variables.emplace_back(0);

    for (auto & coef : coefs) {
        coef.push_back(-coef[idx]);
    }
}

LpProblem::LpProblem(const ObjectiveType& obj, const unsigned int &num_vars) {
    objective = obj;
    obj_coefs = std::vector<double>(num_vars);
    variables = std::vector<LpVariable>(num_vars);
    obj_b = 0;

    num_variables = num_vars;
    objective_inited = false;
}

LpProblem::LpProblem(const ObjectiveType& obj, const Expression &objective_func) {
    objective = obj;
    num_variables = objective_func.arr.size();
    obj_coefs = std::vector<double>(num_variables);
    variables = std::vector<LpVariable>(num_variables);

    for (int i = 0; i < num_variables; ++i)
        obj_coefs[i] = objective_func.arr[i];
    obj_b = objective_func.b;
    objective_inited = true;
}

void LpProblem::add_constraint(const Expression &constraint) {
    if (constraint.exp_type == Expression::NO_SYMBOL) {
        std::cerr << "Expression does not have a valid symbol and will be viewed as a objective." << std::endl;
        set_objective(constraint);
    } else {
        _add_constraint(constraint);
    }
}

void LpProblem::set_objective(const Expression &objective_func) {
    for (int i = 0; i < num_variables; ++i)
        obj_coefs[i] = objective_func.arr[i];
    obj_b = objective_func.b;
    objective_inited = true;
}

const std::vector<LpVariable>& LpProblem::get_variables() const {
    return variables;
}

void LpProblem::set_variable_bound(const int &id, const double &low, const double &up) {
    variables[id].set_bound(low, up);
}

void LpProblem::set_variable_lower_bound(const int &id, const double &low) {
    variables[id].set_lower_bound(low);
}

void LpProblem::set_variable_upper_bound(const int &id, const double &up) {
    variables[id].set_upper_bound(up);
}

LpProblem& LpProblem::operator+=(const Expression& constraint) {
    if (constraint.exp_type == Expression::NO_SYMBOL) {
        if (objective_inited)
            std::cerr << "Expression does not have a valid symbol and will be viewed as a objective." << std::endl;
        set_objective(constraint);
    } else
        _add_constraint(constraint);

    return *this;
}

void LpProblem::display() const {
    if (!objective_inited)
        throw std::invalid_argument("Objective is not initialized");

    bool first_time = true;

    std::cout << "Objective: " << std::endl << (objective == MAXIMIZE ? "MAXIMIZE " : "MINIMIZE ");
    for (int i = 0; i < num_variables; ++i) {
        if (obj_coefs[i] == 0)
            continue;

        if (first_time) {
            first_time = false;
            std::cout << obj_coefs[i] << " * x" << i;
        } else {
            if (obj_coefs[i] > 0)
                std::cout << " + " << obj_coefs[i] << " * x" << i;
            else
                std::cout << " - " << -obj_coefs[i] << " * x" << i;
        }
    }
    if (obj_b > 0)
        std::cout << " + " << obj_b;
    else if (obj_b < 0)
        std::cout << " - " << -obj_b;
    std::cout << std::endl << std::endl;

    std::cout << "Subject to: " << std::endl;
    for (int row = 0; row < coefs.size(); ++row) {
        first_time = true;

        for (int col = 0; col < num_variables; ++col) {
            if (coefs[row][col] == 0)
                continue;

            if (first_time) {
                first_time = false;
                std::cout << coefs[row][col] << " * x" << col;
            } else {
                if (coefs[row][col] > 0)
                    std::cout << " + " << coefs[row][col] << " * x" << col;
                else
                    std::cout << " - " << -coefs[row][col] << " * x" << col;
            }
        }
        std::cout << " <= " << b[row] << std::endl;
    }

    for (int i = 0; i < num_variables; ++i) {
        const bool has_low = variables[i].lower_bound != -std::numeric_limits<double>::infinity();

        if (const bool has_up = variables[i].upper_bound != std::numeric_limits<double>::infinity(); has_low && has_up)
            std::cout << variables[i].lower_bound << " <= x" << i << " <= " << variables[i].upper_bound << std::endl;
        else if (has_low)
            std::cout << 'x' << i  << " >= " << variables[i].lower_bound << std::endl;
        else if (has_up)
            std::cout << "x" << i << variables[i].upper_bound << std::endl;
    }
    std::cout << std::endl;
}
