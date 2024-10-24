# Optimization Models

This repository contains implementations of optimization models, primarily focusing on linear programming using the Simplex Method. Future enhancements will include additional solving methods, such as the Interior Point Method, as well as other optimization models.

## Directory Structure

- **LPtest/**: Contains examples of using the LP model.
- **test/**: Contains test files for various classes.
- **CMakeLists.txt**: Configuration file for CMake.
- **LpProblem.cpp / .h**: Implementation of the linear programming problem class.
- **LpVariable.cpp / .h**: Implementation of the linear programming variable class.
- **Matrix.cpp / .h**: Implementation of matrix operations.
- **SimplexAlgorithm.cpp / .h**: Implementation of the Simplex Algorithm.
- **mainwindow.cpp / .h / .ui**: Implementation of the main window's functionality and UI design.
- **main.cpp**: The main file for the GUI interface.

## Getting Started

To compile the project, use CMake:

```bash
mkdir build
cd build
cmake ..
make
```

## Calling the Model
This optimization model supports two methods of invocation:

### 1. Including in C++ Code
You can directly include the model in your C++ projects. Here is a sample code snippet demonstrating how to use the model (the code snippet is the same as `LP1.cpp`):

```cpp
#include "../SimplexAlgorithm.h"

int main() {
    LpProblem problem(LpProblem::MAXIMIZE, 3);
    problem.set_lower_bounds(0);
    const auto vars = problem.get_variables();

    problem += vars[0] * 6 + vars[1] * 14 + vars[2] * 13;
    problem += vars[0] * 0.5 + vars[1] * 2 + vars[2] <= 24;
    problem += vars[1] * 2 + vars[0] + vars[2] * 4 <= 60;
    problem.set_variable_lower_bound(0, 10);

    problem.display();

    auto simplex = SimplexAlgorithm(problem);
    double optimal_solution = simplex.solve();
    std::cout << "Optimal solution: " << optimal_solution;
}
```

### 2. Using the GUI Interface

We have developed a user-friendly GUI interface that allows users to select the type of problem to solve and the method to use. Users can also choose a CSV file for input. An example CSV file is provided in the `test/` directory as `LpProblem.csv`.

- Select the problem type (e.g., linear programming).
- Choose the solving method (e.g., Simplex Method).
- Load your CSV file through the file explorer.
- This GUI aims to simplify the process of solving optimization problems for users who prefer a visual interface.

## Future Work

- Add more solving methods, including the Interior Point Method.
- Implement additional optimization models.

## Contributing

Feel free to open issues or submit pull requests for improvements or additional features.
