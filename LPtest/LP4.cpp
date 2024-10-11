#include "../SimplexAlgorithm.h"

int main() {
    LpProblem problem(LpProblem::MINIMIZE, 5);
    const auto x = problem.get_variables();

    problem += 40 * x[0] + 35 * x[1] + 15 * x[2] + 18 * x[3] + 10 * x[4];
    problem += 54 * x[0] + 20 * x[1] + 15 * x[2] + 8 * x[3] + 40 * x[4] >= 280;
    problem += 54 * x[0] + 20 * x[1] + 15 * x[2] + 8 * x[3] + 40 * x[4] <= 320;
    problem += 19 * x[0] + 15 * x[1] + 10 * x[2] <= (54 * x[0] + 20 * x[1] + 15 * x[2] + 8 * x[3] + 40 * x[4]) * 0.3;
    problem += 15 * x[2] + 350 * x[3] >= 600;
    problem += x[1] + 3 * x[2] + x[3] >= 10;
    problem += 8 * x[0] + x[2] + x[3] + x[4] >= 30;
    problem += x[1] - 0.5 * x[0] >= 0;
    problem.set_variable_lower_bound(0, 2);

    problem.display();

    auto simplex = SimplexAlgorithm(problem);
    simplex.solve();
}