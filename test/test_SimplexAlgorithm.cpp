#include "../SimplexAlgorithm.h"

int main() {
    std::cout << "The original LP problem: " << std::endl;
    LpProblem problem(LpProblem::MAXIMIZE, 2);

    const auto vars = problem.get_variables();

    problem.set_variable_upper_bound(0, 4);
    problem.set_variable_upper_bound(1, 6);

    const Expression objective = vars[0] * 300 + vars[1] * 500;
    const Expression constraint = vars[0] * 3 + vars[1] * 2 <= 18;

    problem += objective;
    problem += constraint;

    problem.display();

    auto simplex = SimplexAlgorithm(problem);
    simplex.solve();
}