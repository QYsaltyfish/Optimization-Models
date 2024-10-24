//
// Created by QY.
//

#include "../SimplexAlgorithm.h"

int main() {
    LpProblem problem(LpProblem::MAXIMIZE, 2);
    problem.set_lower_bounds(0);
    const auto vars = problem.get_variables();
    const auto& x1 = vars[0];
    const auto& x2 = vars[1];

    problem += x1 + x2;
    problem += x1 + x2 == 10;
    problem += 2 * x1 + x2 >= 2;
    problem.set_variable_lower_bound(0, -std::numeric_limits<double>::infinity());

    problem.display();

    auto simplex = SimplexAlgorithm(problem);
    double optimal_solution = simplex.solve();
    std::cout << "Optimal solution: " << optimal_solution;
}