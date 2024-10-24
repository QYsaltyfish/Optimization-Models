//
// Created by QY.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cctype>
#include "LpProblem.h"

// Implementation of adding a constraint
void LpProblem::_add_constraint(const Expression &constraint) {
    int count = 0;  // Count of non-zero coefficients
    int id = -1;  // Index of the last non-zero coefficient

    // Loop through the number of variables to find non-zero coefficients
    for (int i = 0; i < num_variables; ++i) {
        if (constraint.arr[i] == 0)
            continue;
        ++count;
        id = i;
        if (count == 2)
            break;
    }

    // If two non-zero coefficients are found
    if (count == 2) {
        coefs.push_back(constraint.arr);  // Add coefficients to the matrix
        b.push_back(constraint.b);  // Add right-hand side value

        if (constraint.exp_type == Expression::EQUAL) {
            auto new_vector = std::vector<double>(num_variables);
            for (int i = 0; i < num_variables; ++i)
                new_vector[i] = -constraint.arr[i];  // Negate the coefficients
            coefs.push_back(new_vector);  // Add negated coefficients
            b.push_back(-constraint.b);  // Negate the right-hand side value
        }
    } else {
        if (id == -1)
            return;  // No non-zero coefficients found, exit early
        // Adjust variable bounds based on the constraint
        if (constraint.arr[id] > 0 && constraint.arr[id] < variables[id].upper_bound)
            set_variable_upper_bound(id, constraint.b / constraint.arr[id]);
        else if (constraint.arr[id] < 0 && constraint.arr[id] > variables[id].lower_bound)
            set_variable_lower_bound(id, constraint.b / constraint.arr[id]);
    }
}

// Method to transform the lower bound of a variable
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
            b[row] -= coefs[row][idx] * variables[idx].lower_bound;  // Adjust right-hand side
        }
        obj_b += obj_coefs[idx] * variables[idx].lower_bound;  // Adjust objective constant
        variables[idx].lower_bound = 0;  // Set lower bound to 0
        variables[idx].upper_bound -= variables[idx].lower_bound;  // Adjust upper bound
        return;
    }

    // Case 2: The variable doesn't have a lower bound
    // Case 2.1: The variable has an upper bound
    if (variables[idx].upper_bound != inf) {
        for (int row = 0; row < coefs.size(); ++row) {
            if (coefs[row][idx] == 0)
                continue;
            b[row] -= coefs[row][idx] * variables[idx].upper_bound;  // Adjust right-hand side
            coefs[row][idx] = -coefs[row][idx];  // Negate coefficients
        }
        obj_b += obj_coefs[idx] * variables[idx].upper_bound;  // Adjust objective constant
        variables[idx].lower_bound = 0;
        variables[idx].upper_bound = inf;
        return;
    }

    // Case 2.2: The variable doesn't have an upper bound
    ++num_variables;

    obj_coefs.push_back(0);  // Add a new coefficient for the new variable
    variables.emplace_back(0);  // Add a new variable

    // Adjust coefficients for existing constraints
    for (auto & coef : coefs) {
        coef.push_back(-coef[idx]);
    }
}

// Constructor to initialize from an objective type and number of variables
LpProblem::LpProblem(const ObjectiveType& obj, const unsigned int &num_vars) {
    LpVariable::id_counter = 0;  // Reset variable ID counter
    objective = obj;
    obj_coefs = std::vector<double>(num_vars);
    variables = std::vector<LpVariable>(num_vars);
    obj_b = 0;

    num_variables = num_vars;
    objective_inited = false;
}

// Constructor to initialize from an objective type and objective function
LpProblem::LpProblem(const ObjectiveType& obj, const Expression &objective_func) {
    LpVariable::id_counter = 0;  // Reset variable ID counter
    objective = obj;
    num_variables = objective_func.arr.size();
    obj_coefs = std::vector<double>(num_variables);
    variables = std::vector<LpVariable>(num_variables);

    for (int i = 0; i < num_variables; ++i)
        obj_coefs[i] = objective_func.arr[i];
    obj_b = objective_func.b;
    objective_inited = true;
}

// Constructor to initialize from a CSV file
LpProblem::LpProblem(const std::string &csvFilePath) {
    LpVariable::id_counter = 0;  // Reset variable ID counter
    std::ifstream file(csvFilePath);  // Open the CSV file
    std::string line;  // Line buffer for reading

    // Read the first line to determine the objective type
    if (std::getline(file, line)) {
        std::transform(line.begin(), line.end(), line.begin(),
                       [](unsigned char c) { return std::toupper(c); });  // Convert to uppercase
        if (line == "MAX" || line == "MAXIMIZE")
            objective = MAXIMIZE;
        else if (line == "MIN" || line == "MINIMIZE")
            objective = MINIMIZE;
        else
            throw std::invalid_argument("The objective is not supported: " + line);  // Invalid objective
    } else {
        throw std::invalid_argument("The first line must be the type of objective.");  // Empty file
    }

    num_variables = 0;
    obj_coefs = std::vector<double>();

    // Read the second line for the objective function coefficients
    if (std::getline(file, line)) {
        std::stringstream ss(line);  // String stream for parsing
        std::string token;

        // Parse coefficients
        while (std::getline(ss, token, ',')) {
            ++num_variables;
            obj_coefs.push_back(std::stod(token));  // Add coefficient
        }
        --num_variables;  // Adjust variable count for constant
        obj_b = obj_coefs[num_variables];  // Set objective constant
        obj_coefs.pop_back();  // Remove the constant from coefficients

        std::cout << std::endl;
    } else {
        throw std::invalid_argument("The second line must be the objective function.");  // Missing objective function
    }
    variables = std::vector<LpVariable>(num_variables);  // Initialize variables

    enum flagType {UNKNOWN, LESS_EQUAL, EQUAL, GREATER_EQUAL};  // Enumeration for constraint types

    // Read constraints from the CSV file
    while (std::getline(file, line)) {
        std::stringstream ss(line);  // String stream for parsing
        std::string token;
        flagType flag;
        bool b_inited = false;  // Flag for initializing 'b'
        int col = 0;
        int row_coef_count = 0;  // Count of non-zero coefficients in the current row
        int row_coef_idx;  // Index of the last non-zero coefficient

        flag = UNKNOWN;
        coefs.emplace_back(num_variables);  // Create a new row for coefficients
        while (std::getline(ss, token, ',')) {
            if (flag != UNKNOWN) {
                if (flag == LESS_EQUAL) {
                    b.push_back(std::stod(token));  // Add right-hand side for <=
                } else if (flag == GREATER_EQUAL) {
                    b.push_back(-std::stod(token));  // Add negated right-hand side for >=
                    for (double &coef: coefs[coefs.size() - 1])
                        coef = -coef;  // Negate coefficients
                } else {
                    b.push_back(std::stod(token));  // Add right-hand side for =
                    coefs.emplace_back(num_variables);  // Create a new row for coefficients
                    for (int i = 0; i < num_variables; ++i) {
                        coefs[coefs.size() - 1][i] = -coefs[coefs.size() - 2][i];  // Negate previous coefficients
                    }
                    b.push_back(-std::stod(token));  // Add negated right-hand side
                }

                b_inited = true;
                continue;
            }

            // Determine the constraint type from the token
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
                throw std::invalid_argument("Too many constraint coefficients.");  // Invalid number of coefficients

            double coef = std::stod(token);
            if (coef != 0) {
                ++row_coef_count;
                row_coef_idx = col;
            }
            coefs[coefs.size() - 1][col] = coef;  // Add coefficient to the matrix
            ++col;
        }

        if (flag == UNKNOWN) {
            throw std::invalid_argument("No <=/==/>= in at least one line.");  // No constraint type found
        }
        if (!b_inited) {
            throw std::invalid_argument("b is not initialized in at least one line.");  // 'b' not initialized
        }

        if (row_coef_count == 0) {
            throw std::invalid_argument("All coefs are 0 in at least one line.");  // No non-zero coefficients
        }

        // Handle the case when there is only one non-zero coefficient
        if (row_coef_count == 1) {
            if (flag != EQUAL) {
                double coef = coefs[coefs.size() - 1][row_coef_idx];
                if (coef > 0) {
                    variables[row_coef_idx].set_upper_bound(b[coefs.size() - 1] / coef);
                } else {
                    variables[row_coef_idx].set_lower_bound(b[coefs.size() - 1] / coef);
                }

                b.pop_back();  // Remove the last element from 'b'
                coefs.pop_back();  // Remove the last row from coefficients
            } else {
                // TODO: Handle when flag is EQUAL
                throw std::invalid_argument("The situation when a variable equals a constant is not supported yet.");
            }
        }
    }

    objective_inited = true;
}

// Method to add a constraint or set it as the objective
void LpProblem::add_constraint(const Expression &constraint) {
    LpVariable::id_counter = 0;  // Reset variable ID counter
    if (constraint.exp_type == Expression::NO_SYMBOL) {
        std::cerr << "Expression does not have a valid symbol and will be viewed as an objective." << std::endl;
        set_objective(constraint);
    } else {
        _add_constraint(constraint);
    }
}

// Method to set the objective function
void LpProblem::set_objective(const Expression &objective_func) {
    for (int i = 0; i < num_variables; ++i)
        obj_coefs[i] = objective_func.arr[i];
    obj_b = objective_func.b;
    objective_inited = true;
}

// Getter method for variables
const std::vector<LpVariable>& LpProblem::get_variables() const {
    return variables;
}

// Set bounds for a specific variable
void LpProblem::set_variable_bound(const int &id, const double &low, const double &up) {
    variables[id].set_bound(low, up);
}

// Set the same bounds for all variables
void LpProblem::set_bounds(const double &low, const double &up) {
    for (int i = 0; i < num_variables; ++i)
        variables[i].set_bound(low, up);
}

// Set lower bound for a specific variable
void LpProblem::set_variable_lower_bound(const int &id, const double &low) {
    variables[id].set_lower_bound(low);
}

// Set the same lower bound for all variables
void LpProblem::set_lower_bounds(const double &low) {
    for (int i = 0; i < num_variables; ++i) {
        variables[i].set_lower_bound(low);
    }
}

// Set upper bound for a specific variable
void LpProblem::set_variable_upper_bound(const int &id, const double &up) {
    variables[id].set_upper_bound(up);
}

// Set the same upper bound for all variables
void LpProblem::set_upper_bounds(const double &up) {
    for (int i = 0; i < num_variables; ++i) {
        variables[i].set_upper_bound(up);
    }
}

// Overloaded operator to add constraints
LpProblem& LpProblem::operator+=(const Expression& constraint) {
    if (constraint.exp_type == Expression::NO_SYMBOL) {
        if (objective_inited)
            std::cerr << "Expression does not have a valid symbol and will be viewed as an objective." << std::endl;
        set_objective(constraint);
    } else
        _add_constraint(constraint);

    return *this;
}

// Display the linear programming problem in a readable format
void LpProblem::display() const {
    if (!objective_inited)
        throw std::invalid_argument("Objective is not initialized");  // Objective not set

    bool first_time = true;  // Flag to handle the first term in the objective function

    std::cout << "Objective: " << std::endl << (objective == MAXIMIZE ? "MAXIMIZE " : "MINIMIZE ");
    for (int i = 0; i < num_variables; ++i) {
        if (obj_coefs[i] == 0)
            continue;  // Skip zero coefficients

        if (first_time) {
            first_time = false;
            std::cout << obj_coefs[i] << " * x" << i;  // Print first term
        } else {
            if (obj_coefs[i] > 0)
                std::cout << " + " << obj_coefs[i] << " * x" << i;  // Print positive term
            else
                std::cout << " - " << -obj_coefs[i] << " * x" << i;  // Print negative term
        }
    }
    if (obj_b > 0)
        std::cout << " + " << obj_b;  // Add constant if positive
    else if (obj_b < 0)
        std::cout << " - " << -obj_b;  // Add constant if negative
    std::cout << std::endl << std::endl;

    std::cout << "Subject to: " << std::endl;  // Print constraints
    for (int row = 0; row < coefs.size(); ++row) {
        first_time = true;

        for (int col = 0; col < num_variables; ++col) {
            if (coefs[row][col] == 0)
                continue;

            if (first_time) {
                first_time = false;
                std::cout << coefs[row][col] << " * x" << col;  // Print first term
            } else {
                if (coefs[row][col] > 0)
                    std::cout << " + " << coefs[row][col] << " * x" << col;  // Print positive term
                else
                    std::cout << " - " << -coefs[row][col] << " * x" << col;  // Print negative term
            }
        }
        std::cout << " <= " << b[row] << std::endl;  // Print right-hand side
    }

    // Print variable bounds
    for (int i = 0; i < num_variables; ++i) {
        const bool has_low = variables[i].lower_bound != -inf;
        const bool has_up = variables[i].upper_bound != inf;

        if (has_low && has_up)
            // Both bounds
            std::cout << variables[i].lower_bound << " <= x" << i << " <= " << variables[i].upper_bound << std::endl;
        else if (has_low)
            std::cout << 'x' << i << " >= " << variables[i].lower_bound << std::endl;  // Lower bound only
        else if (has_up)
            std::cout << 'x' << i << " <= " << variables[i].upper_bound << std::endl;  // Upper bound only
    }
    std::cout << std::endl;  // Newline for readability
}
