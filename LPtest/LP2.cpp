#include "../SimplexAlgorithm.h"

int main() {
    LpProblem problem(LpProblem::MINIMIZE, 2);
    const auto vars = problem.get_variables();

    problem += 3 * vars[0] + 2 * vars[1];
    problem += vars[0] + 2 * vars[1] <= 12;
    problem += 2 * vars[0] + 3 * vars[1] == 12;
    problem += 2 * vars[0] + vars[1] >= 8;

    problem.display();

    auto simplex = SimplexAlgorithm(problem);
    simplex.solve();
}