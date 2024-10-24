#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cctype>
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
    if (variables[idx].lower_bound != -inf) {

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
    if (variables[idx].upper_bound != inf) {
        for (int row = 0; row < coefs.size(); ++row) {
            if (coefs[row][idx] == 0)
                continue;
            b[row] -= coefs[row][idx] * variables[idx].upper_bound;
            coefs[row][idx] = -coefs[row][idx];
        }
        obj_b += obj_coefs[idx] * variables[idx].upper_bound;
        variables[idx].lower_bound = 0;
        variables[idx].upper_bound = inf;
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
    LpVariable::id_counter = 0;
    objective = obj;
    obj_coefs = std::vector<double>(num_vars);
    variables = std::vector<LpVariable>(num_vars);
    obj_b = 0;

    num_variables = num_vars;
    objective_inited = false;
}

LpProblem::LpProblem(const ObjectiveType& obj, const Expression &objective_func) {
    LpVariable::id_counter = 0;
    objective = obj;
    num_variables = objective_func.arr.size();
    obj_coefs = std::vector<double>(num_variables);
    variables = std::vector<LpVariable>(num_variables);

    for (int i = 0; i < num_variables; ++i)
        obj_coefs[i] = objective_func.arr[i];
    obj_b = objective_func.b;
    objective_inited = true;
}

LpProblem::LpProblem(const std::string &csvFilePath) {
    LpVariable::id_counter = 0;
    std::ifstream file(csvFilePath);
    std::string line;

    if (std::getline(file, line)) {
        std::transform(line.begin(), line.end(), line.begin(),
                       [](unsigned char c) { return std::toupper(c); });
        if (line == "MAX" || line == "MAXIMIZE")
            objective = MAXIMIZE;
        else if (line == "MIN" || line == "MINIMIZE")
            objective = MINIMIZE;
        else
            throw std::invalid_argument("The objective is not supported: " + line);
    } else {
        throw std::invalid_argument("The first line must be the type of objective.");
    }

    num_variables = 0;
    obj_coefs = std::vector<double>();
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;

        while (std::getline(ss, token, ',')) {
            ++num_variables;
            obj_coefs.push_back(std::stod(token));
        }
        --num_variables;
        obj_b = obj_coefs[num_variables];
        obj_coefs.pop_back();

        std::cout << std::endl;
    } else {
        throw std::invalid_argument("The second line must be the objective function.");
    }
    variables = std::vector<LpVariable>(num_variables);

    enum flagType {UNKNOWN, LESS_EQUAL, EQUAL, GREATER_EQUAL};
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        flagType flag;
        bool b_inited = false;
        int col = 0;
        int row_coef_count = 0;
        int row_coef_idx;

        flag = UNKNOWN;
        coefs.emplace_back(num_variables);
        while (std::getline(ss, token, ',')) {
            if (flag != UNKNOWN) {
                if (flag == LESS_EQUAL) {
                    b.push_back(std::stod(token));
                } else if (flag == GREATER_EQUAL) {
                    b.push_back(-std::stod(token));
                    for (double &coef: coefs[coefs.size() - 1])
                        coef = -coef;
                } else {
                    b.push_back(std::stod(token));
                    coefs.emplace_back(num_variables);
                    for (int i = 0; i < num_variables; ++i) {
                        coefs[coefs.size() - 1][i] = -coefs[coefs.size() - 2][i];
                    }
                    b.push_back(-std::stod(token));
                }

                b_inited = true;
                continue;
            }

            if (token == "<" || token == "<=") {
                flag = LESS_EQUAL;
                continue;
            } else if (token == "=" || token == "==") {
                flag = EQUAL;
                continue;
            } else if (token == ">" || token == ">=") {
                flag = GREATER_EQUAL;
                continue;
            }

            if (col >= num_variables)
                throw std::invalid_argument("Too many constraint coefficients.");

            double coef = std::stod(token);
            if (coef != 0) {
                ++row_coef_count;
                row_coef_idx = col;
            }
            coefs[coefs.size() - 1][col] = coef;
            ++col;
        }

        if (flag == UNKNOWN) {
            throw std::invalid_argument("No <=/==/>= in at least one line.");
        }
        if (!b_inited) {
            throw std::invalid_argument("b is not initialized in at least one line.");
        }

        if (row_coef_count == 0) {
            throw std::invalid_argument("All coefs are 0 in at least one line.");
        }

        if (row_coef_count == 1) {
            if (flag != EQUAL) {
                double coef = coefs[coefs.size() - 1][row_coef_idx];
                if (coef > 0) {
                    variables[row_coef_idx].set_upper_bound(b[coefs.size() - 1] / coef);
                } else {
                    variables[row_coef_idx].set_lower_bound(b[coefs.size() - 1] / coef);
                }

                b.pop_back();
                coefs.pop_back();
            } else {
                // TODO: Handle when flag is EQUAL
                throw std::invalid_argument("The situation when a variable equals a constant is not supported yet.");
            }
        }
    }

    objective_inited = true;
}

void LpProblem::add_constraint(const Expression &constraint) {
    LpVariable::id_counter = 0;
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

void LpProblem::set_bounds(const double &low, const double &up) {
    for (int i = 0; i < num_variables; ++i)
        variables[i].set_bound(low, up);
}

void LpProblem::set_variable_lower_bound(const int &id, const double &low) {
    variables[id].set_lower_bound(low);
}

void LpProblem::set_lower_bounds(const double &low) {
    for (int i = 0; i < num_variables; ++i) {
        variables[i].set_lower_bound(low);
    }
}

void LpProblem::set_variable_upper_bound(const int &id, const double &up) {
    variables[id].set_upper_bound(up);
}

void LpProblem::set_upper_bounds(const double &up) {
    for (int i = 0; i < num_variables; ++i) {
        variables[i].set_upper_bound(up);
    }
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
        const bool has_low = variables[i].lower_bound != -inf;
        const bool has_up = variables[i].upper_bound != inf;

        if (has_low && has_up)
            std::cout << variables[i].lower_bound << " <= x" << i << " <= " << variables[i].upper_bound << std::endl;
        else if (has_low)
            std::cout << 'x' << i << " >= " << variables[i].lower_bound << std::endl;
        else if (has_up)
            std::cout << 'x' << i << " <= " << variables[i].upper_bound << std::endl;
    }
    std::cout << std::endl;
}
