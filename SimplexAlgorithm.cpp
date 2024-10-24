//
// Created by QY.
//

#include "SimplexAlgorithm.h"

// SimplexAlgorithm class implementation for solving linear programming problems

// Check for the entering column in the tableau
int SimplexAlgorithm::check_column() {
    // Create a vector to hold the coefficients of the basic variables
    auto basis_coefs = std::vector<double>(lp.coefs.size());

    // Fill basis_coefs with the coefficients of the basic variables
    for (int i = 0; i < basis_coefs.size(); ++i) {
        basis_coefs[i] = lp.obj_coefs[basis[i]];
    }

    double max = 0;
    int max_idx = -1;

    // Loop through all variables to find the entering variable
    for (int i = 0; i < lp.num_variables; ++i) {
        if (is_basis[i]) {
            continue;  // Skip basic variables
        }

        double sigma = 0;  // Calculate the net benefit of entering variable
        for (int j = 0; j < lp.coefs.size(); ++j) {
            sigma += basis_coefs[j] * lp.coefs[j][i];
        }

        sigma = lp.obj_coefs[i] - sigma;  // Calculate the reduced cost

        // Update maximum reduced cost and index if necessary
        if (sigma > max) {
            max = sigma;
            max_idx = i;
        }
    }

    return max_idx;  // Return the index of the entering variable
}

// Check for the leaving row in the tableau based on the entering column index
int SimplexAlgorithm::check_row(const int& idx) {
    double min = inf;
    int min_idx = -1;

    // Loop through each constraint to find the leaving variable
    for (int row = 0; row < lp.coefs.size(); ++row) {
        if (lp.coefs[row][idx] > 0) { // Only consider positive coefficients
            // Calculate the ratio of right-hand side to coefficient
            if (const double curr_ratio = lp.b[row] / lp.coefs[row][idx]; curr_ratio < min) {
                min = curr_ratio;
                min_idx = row;
            }
        }
    }
    return min_idx;  // Return the index of the leaving variable
}

// Solve the linear programming problem using the Simplex method
double SimplexAlgorithm::solve() {
    int state;  // State of the solution process

    while (true) {
        int idx = check_column();  // Check for entering column

        // If no entering column is found, the problem is optimal
        if (idx == -1) {
            state = 0;
            break;
        }

        const int row = check_row(idx);  // Check for leaving row

        // If no leaving row is found, the problem is unbounded
        if (row == -1) {
            state = -1;
            break;
        }

        double coef_temp = lp.coefs[row][idx];

        // Normalize the pivot row
        for (auto& coef: lp.coefs[row])
            coef /= coef_temp;
        lp.b[row] /= coef_temp;

        // Update other rows
        for (int i = 0; i < row; ++i) {
            coef_temp = lp.coefs[i][idx];
            for (int col = 0; col < lp.num_variables; ++col) {
                lp.coefs[i][col] -= coef_temp * lp.coefs[row][col];  // Row operation
            }
            lp.b[i] -= coef_temp * lp.b[row];  // Update right-hand side
        }
        for (int i = row + 1; i < lp.coefs.size(); ++i) {
            coef_temp = lp.coefs[i][idx];
            for (int col = 0; col < lp.num_variables; ++col) {
                lp.coefs[i][col] -= coef_temp * lp.coefs[row][col];  // Row operation
            }
            lp.b[i] -= coef_temp * lp.b[row];  // Update right-hand side
        }

        // Update basis
        is_basis[basis[row]] = false;
        basis[row] = idx;
        is_basis[idx] = true;
    }

    // Calculate the optimal value
    double optimal_value = 0;
    for (int i = 0; i < lp.coefs.size(); ++i) {
        optimal_value += lp.b[i] * lp.obj_coefs[basis[i]];
    }
    optimal_value += lp.obj_b;
    optimal_value = (reversed ? -optimal_value : optimal_value);

    // Store the solution in the solution vector
    for (int i = 0; i < lp.b.size(); ++i) {
        solution[basis[i]] = lp.b[i];
    }

    // Return the state of the solution
    if (state == -1)
        return inf;  // Unbounded case
    else {
        return optimal_value;  // Optimal value
    }
}

// Constructor for the SimplexAlgorithm class
SimplexAlgorithm::SimplexAlgorithm(const LpProblem &problem): lp(problem), reversed(false) {
    basis = std::vector<unsigned int>(lp.coefs.size());
    is_basis = std::vector<bool>(lp.num_variables, false);
    transform_problem();  // Transform the problem for Simplex
    solution = std::vector<double>(lp.num_variables, 0);
    lp.display();
}

// Transform the problem to standard form
void SimplexAlgorithm::transform_problem() {
    for (int i = 0; i < lp.num_variables; ++i) {
        lp.transform_lower_bound(i);
        transform_upper_bound(i);
    }

    for (int i = 0; i < lp.coefs.size(); ++i) {
        transform_less_equal(i);  // Transform inequalities
    }

    // Adjust for negative right-hand side values
    for (int i = 0; i < lp.coefs.size(); ++i) {
        if (lp.b[i] < 0) {
            lp.b[i] = -lp.b[i];
            for (double & col : lp.coefs[i])
                col = -col;
        }
    }

    // If the objective is minimization, convert to maximization
    if (lp.objective == LpProblem::MINIMIZE) {
        lp.objective = LpProblem::MAXIMIZE;
        reversed = true;
        for (int i = 0; i < lp.num_variables; ++i)
            lp.obj_coefs[i] = -lp.obj_coefs[i];  // Negate objective coefficients
        lp.obj_b = -lp.obj_b;  // Negate the objective value
    }
}

// Transform LESS_EQUAL constraints to EQUAL by adding slack variables
void SimplexAlgorithm::transform_less_equal(const unsigned int &idx) {
    if (lp.b[idx] >= 0) {
        // Add slack variable for non-negative constraint
        for (unsigned int row = 0; row < idx; ++row)
            lp.coefs[row].push_back(0);  // Add zero column for slack variable
        for (unsigned int row = idx + 1; row < lp.coefs.size(); ++row)
            lp.coefs[row].push_back(0);  // Add zero column for slack variable
        lp.coefs[idx].push_back(1);  // Add slack variable
        lp.variables.emplace_back(0);  // Add slack variable to variable list
        lp.obj_coefs.push_back(0);  // Objective coefficient for slack variable

        // Update basis
        basis[idx] = lp.num_variables;
        is_basis[lp.num_variables] = true;
        ++lp.num_variables;
    } else {
        // Introduce artificial variables for negative constraints
        for (auto& coef: lp.coefs[idx])
            coef = -coef;  // Negate coefficients
        lp.b[idx] = -lp.b[idx];  // Negate right-hand side

        // Adjust tableau for artificial variables
        for (unsigned int row = 0; row < idx; ++row) {
            lp.coefs[row].push_back(0);
            lp.coefs[row].push_back(0);
        }
        for (unsigned int row = idx + 1; row < lp.coefs.size(); ++row) {
            lp.coefs[row].push_back(0);
            lp.coefs[row].push_back(0);
        }
        lp.coefs[idx].push_back(-1);  // Add slack variable
        lp.coefs[idx].push_back(1);  // Add artificial variable
        lp.variables.emplace_back(0);  // Add slack variable to variable list
        lp.variables.emplace_back(0);  // Add artificial variable to variable list
        lp.obj_coefs.push_back(0);  // Objective coefficient for slack variable

        // Add a large constant for artificial variable in the objective
        if (lp.objective == LpProblem::MAXIMIZE)
            lp.obj_coefs.push_back(-10000000);
        else
            lp.obj_coefs.push_back(10000000);

        // Update basis
        basis[idx] = lp.num_variables + 1;  // Set basis for artificial variable
        is_basis[lp.num_variables + 1] = true;

        lp.num_variables += 2;
    }

}

// Transform upper bounds by adding a surplus variable
void SimplexAlgorithm::transform_upper_bound(const int &idx) {
    // If there is an upper bound, we need to adjust the model
    if (lp.variables[idx].upper_bound != inf) {
        lp.variables.emplace_back(0);  // Add a surplus variable to the variable list
        basis.push_back(lp.num_variables);  // Update basis
        ++lp.num_variables;

        // Adjust coefficients in the tableau
        for (auto & coef : lp.coefs) {
            coef.push_back(0);  // Add a column for surplus variable
        }

        // Create a new vector for the new constraint row
        auto new_vector = std::vector<double>(lp.num_variables, 0);
        new_vector[lp.num_variables - 1] = 1;  // Coefficient for surplus variable
        new_vector[idx] = 1;  // Coefficient for the original variable
        lp.coefs.push_back(new_vector);

        lp.b.push_back(lp.variables[idx].upper_bound);  // Update right-hand side
        lp.variables[idx].upper_bound = inf;
        lp.obj_coefs.push_back(0);  // Objective coefficient for surplus variable

        is_basis.push_back(true);  // Mark the surplus variable as basic
    }
}
