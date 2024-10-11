#include "SimplexAlgorithm.h"


int SimplexAlgorithm::check_column() {
    auto basis_coefs = std::vector<double>(lp.coefs.size());

    for (int i = 0; i < basis_coefs.size(); ++i) {
        basis_coefs[i] = lp.obj_coefs[basis[i]];
    }

    double max = 0;
    int max_idx = -1;

    for (int i = 0; i < lp.num_variables; ++i) {
        if (is_basis[i]) {
            continue;
        }

        double sigma = 0;
        for (int j = 0; j < lp.coefs.size(); ++j) {
            sigma += basis_coefs[j] * lp.coefs[j][i];
        }

        sigma = lp.obj_coefs[i] - sigma;

        if (sigma > max) {
            max = sigma;
            max_idx = i;
        }
    }

    return max_idx;
}

int SimplexAlgorithm::check_row(const int& idx) {
    double min = std::numeric_limits<double>::infinity();
    int min_idx = -1;

    for (int row = 0; row < lp.coefs.size(); ++row) {
        if (lp.coefs[row][idx] > 0) {
            if (const double curr_ratio = lp.b[row] / lp.coefs[row][idx]; curr_ratio < min) {
                min = curr_ratio;
                min_idx = row;
            }
        }
    }
    return min_idx;
}

void SimplexAlgorithm::solve() {
    int state;

    while (true) {
        int idx = check_column();

        // The problem is optimal as c <= 0
        if (idx == -1) {
            state = 0;
            break;
        }

        const int row = check_row(idx);

        // The problem is unbounded
        if (row == -1) {
            state = -1;
            break;
        }

        double coef_temp = lp.coefs[row][idx];

        // Replace the basis
        for (auto& coef: lp.coefs[row])
            coef /= coef_temp;
        lp.b[row] /= coef_temp;

        for (int i = 0; i < row; ++i) {
            coef_temp = lp.coefs[i][idx];
            for (int col = 0; col < lp.num_variables; ++col) {
                lp.coefs[i][col] -= coef_temp * lp.coefs[row][col];
            }
            lp.b[i] -= coef_temp * lp.b[row];
        }
        for (int i = row + 1; i < lp.coefs.size(); ++i) {
            coef_temp = lp.coefs[i][idx];
            for (int col = 0; col < lp.num_variables; ++col) {
                lp.coefs[i][col] -= coef_temp * lp.coefs[row][col];
            }
            lp.b[i] -= coef_temp * lp.b[row];
        }

        is_basis[basis[row]] = false;
        basis[row] = idx;
        is_basis[idx] = true;
    }

    double optimal_value = 0;
    for (int i = 0; i < lp.coefs.size(); ++i) {
        optimal_value += lp.b[i] * lp.obj_coefs[basis[i]];
    }
    optimal_value += lp.obj_b;
    optimal_value = (reversed ? -optimal_value : optimal_value);

    for (int i = 0; i < lp.b.size(); ++i) {
        solution[basis[i]] = lp.b[i];
    }

    if (state == -1)
        std::cout << "The problem is unbounded." << std::endl;
    else {
        std::cout << "The optimal value is: " << optimal_value << std::endl;
        std::cout << "The solution is: ";
        for (int i = 0; i < lp.num_variables; ++i)
            std::cout << 'x' << i << " = " << solution[i] << ' ';
    }

}

SimplexAlgorithm::SimplexAlgorithm(const LpProblem &problem): lp(problem), reversed(false) {
    basis = std::vector<unsigned int>(lp.coefs.size());
    is_basis = std::vector<bool>(lp.num_variables, false);
    transform_problem();
    solution = std::vector<double>(lp.num_variables, 0);
    lp.display();
}

void SimplexAlgorithm::transform_problem() {
    for (int i = 0; i < lp.num_variables; ++i) {
        lp.transform_lower_bound(i);
        transform_upper_bound(i);
    }

    for (int i = 0; i < lp.coefs.size(); ++i) {
        transform_less_equal(i);
    }

    for (int i = 0; i < lp.coefs.size(); ++i) {
        if (lp.b[i] < 0) {
            lp.b[i] = -lp.b[i];
            for (double & col : lp.coefs[i])
                col = -col;
        }
    }

    if (lp.objective == LpProblem::MINIMIZE) {
        lp.objective = LpProblem::MAXIMIZE;
        reversed = true;
        for (int i = 0; i < lp.num_variables; ++i)
            lp.obj_coefs[i] = -lp.obj_coefs[i];
        lp.obj_b = -lp.obj_b;
    }
}

void SimplexAlgorithm::transform_less_equal(const unsigned int &idx) {

    // Transform LESS_EQUAL to EQUAL by adding slack variable.
    if (lp.b[idx] >= 0) {
        for (unsigned int row = 0; row < idx; ++row)
            lp.coefs[row].push_back(0);
        for (unsigned int row = idx + 1; row < lp.coefs.size(); ++row)
            lp.coefs[row].push_back(0);
        lp.coefs[idx].push_back(1);
        lp.variables.emplace_back(0);
        lp.obj_coefs.push_back(0);

        basis[idx] = lp.num_variables;
        is_basis[lp.num_variables] = true;
        ++lp.num_variables;
    } else {
        // Introduce artificial variables
        for (auto& coef: lp.coefs[idx])
            coef = -coef;
        lp.b[idx] = -lp.b[idx];

        for (unsigned int row = 0; row < idx; ++row) {
            lp.coefs[row].push_back(0);
            lp.coefs[row].push_back(0);
        }
        for (unsigned int row = idx + 1; row < lp.coefs.size(); ++row) {
            lp.coefs[row].push_back(0);
            lp.coefs[row].push_back(0);
        }
        lp.coefs[idx].push_back(-1);
        lp.coefs[idx].push_back(1);
        lp.variables.emplace_back(0);
        lp.variables.emplace_back(0);
        lp.obj_coefs.push_back(0);

        // Add a large constant
        if (lp.objective == LpProblem::MAXIMIZE)
            lp.obj_coefs.push_back(-10000000);
        else
            lp.obj_coefs.push_back(10000000);

        basis[idx] = lp.num_variables + 1;
        is_basis[lp.num_variables + 1] = true;

        lp.num_variables += 2;
    }

}

void SimplexAlgorithm::transform_upper_bound(const int &idx) {
    // Transform upper bounds by adding surplus variable.
    if (lp.variables[idx].upper_bound != std::numeric_limits<double>::infinity()) {
        lp.variables.emplace_back(0);
        basis.push_back(lp.num_variables);
        ++lp.num_variables;

        for (auto & coef : lp.coefs) {
            coef.push_back(0);
        }

        auto new_vector = std::vector<double>(lp.num_variables, 0);
        new_vector[lp.num_variables - 1] = 1;
        new_vector[idx] = 1;
        lp.coefs.push_back(new_vector);

        lp.b.push_back(lp.variables[idx].upper_bound);
        lp.variables[idx].upper_bound = std::numeric_limits<double>::infinity();
        lp.obj_coefs.push_back(0);

        is_basis.push_back(true);
    }
}
