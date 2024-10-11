#include "../LpProblem.h"


int main() {
    LpProblem problem(2);

    std::cout << std::endl << "Generate a LP problem with 2 variables: " << std::endl;
    const auto vars = problem.get_variables();
    std::cout << vars[0] << ' ' << vars[1] << std::endl << std::endl;

    std::cout << "Set bounds for variables: " << std::endl;
    problem.set_variable_upper_bound(0, 4);
    std::cout << "0 <= x0 <= 4" << std::endl << std::endl;

    std::cout << "Construct an objective function: " << std::endl;
    const Expression objective = vars[0] * 300 + vars[1] * 500;
    std::cout << "Objective: " << objective << std::endl << std::endl;

    std::cout << "Construct a constraint: " << std::endl;
    const Expression constraint1 = vars[0] * 3 + vars[1] * 2 <= 18 ;
    std::cout << "Constraint1: " << constraint1 << std::endl << std::endl;

    std::cout << "We can also use constraints to set bounds for variables: " << std::endl;
    const Expression constraint2 = vars[1] * 1 <= 6;
    std::cout << "Constraint2: " << constraint2 << std::endl << std::endl;

    std::cout << "Construct a LP problem and display: " << std::endl;
    problem += objective;
    problem += constraint1;
    problem += constraint2;
    problem.display();

    std::cout << std::endl;
    std::cout << "Now we can transform the problem to standard form: " << std::endl;
    problem.transform_problem();
    problem.display();
}
