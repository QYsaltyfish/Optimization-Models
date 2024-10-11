#include "../SimplexAlgorithm.h"

int main() {
    LpProblem problem(LpProblem::MAXIMIZE, 3);
    const auto vars = problem.get_variables();

    problem += vars[0] * 6 + vars[1] * 14 + vars[2] * 13;
    problem += vars[0] * 0.5 + vars[1] * 2 + vars[2] <= 24;
    problem += vars[1] * 2 + vars[0] + vars[2] * 4 <= 60;
    problem.set_variable_lower_bound(0, 10);

    problem.display();

    auto simplex = SimplexAlgorithm(problem);
    simplex.solve();
}