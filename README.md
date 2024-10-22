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

## Getting Started

To compile the project, use CMake:

```bash
mkdir build
cd build
cmake ..
make
```

## Future Work

- Add more solving methods, including the Interior Point Method.
- Implement additional optimization models.

## Contributing

Feel free to open issues or submit pull requests for improvements or additional features.
