//
// Created by QY.
//

#ifndef LpVariable_H
#define LpVariable_H

#include <iostream>
#include <limits>
#include <vector>

#define inf std::numeric_limits<double>::infinity()

class Expression;
class SimplexAlgorithm;

class LpVariable {
    friend class LpProblem;
    friend class Expression;
    friend class SimplexAlgorithm;
    double lower_bound;
    double upper_bound;

    static int id_counter;
    int id;
    double value;

public:
    explicit LpVariable(const double& low = -inf, const double& up = inf);

    Expression operator+(const LpVariable& other) const;
    Expression operator-(const LpVariable& other) const;
    Expression operator+(const Expression& expr) const;
    Expression operator-(const Expression& expr) const;
    Expression operator*(const double& scalar) const;
    Expression operator/(const double& scalar) const;

    void set_bound(const double& low, const double& up);
    void set_lower_bound(const double& low);
    void set_upper_bound(const double& up);

    friend Expression operator*(const double& scalar, const LpVariable& var);
    friend std::ostream& operator<<(std::ostream& os, const LpVariable& var);
};

class Expression {
    friend class LpProblem;
    friend class SimplexAlgorithm;
    enum ExpType {NO_SYMBOL, LESS_EQUAL, EQUAL};
public:
    std::vector<double> arr;
    double b = 0;
    ExpType exp_type = NO_SYMBOL;

    explicit Expression(const unsigned int& arr_size, const bool& init = true);
    Expression(const Expression& other);

    Expression operator+(const LpVariable& var) const;
    Expression operator+(const Expression& other) const;
    Expression operator+(const double& scalar) const;
    Expression& operator+=(const LpVariable& var);
    Expression& operator+=(const Expression& other);
    Expression& operator+=(const double& scalar);

    Expression operator-(const LpVariable& var) const;
    Expression operator-(const Expression& other) const;
    Expression operator-(const double& scalar) const;
    Expression& operator-=(const LpVariable& var);
    Expression& operator-=(const Expression& other);
    Expression& operator-=(const double& scalar);

    Expression operator*(const double& scalar) const;
    Expression operator/(const double& scalar) const;
    Expression& operator*=(const double& scalar);
    Expression& operator/=(const double& scalar);

    Expression operator<=(const Expression& other) const;
    Expression operator<=(const LpVariable& var) const;
    Expression operator<=(const double& scalar) const;

    Expression operator>=(const Expression& other) const;
    Expression operator>=(const LpVariable& var) const;
    Expression operator>=(const double& scalar) const;

    Expression operator==(const Expression& other) const;
    Expression operator==(const LpVariable& var) const;
    Expression operator==(const double& scalar) const;

    friend std::ostream& operator<<(std::ostream& os, const Expression& expr);
};

#endif
