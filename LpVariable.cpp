#include "LpVariable.h"

int LpVariable::id_counter = 0;

Expression::Expression(const unsigned int &arr_size, const bool& init) {
    if (init)
        arr = std::vector<double>(arr_size);
    else
        arr = std::vector<double>(arr_size, 0);
}

Expression::Expression(const Expression& other) {
    arr = other.arr;
    b = other.b;
    exp_type = other.exp_type;
}

Expression Expression::operator+(const LpVariable& var) const{
    if (arr.size() < var.id)
        throw std::invalid_argument("Size of expression doesn't match variables.");

    Expression result(*this);
    ++result.arr[var.id];

    return result;
}

Expression Expression::operator+(const Expression& other) const{
    if (arr.size() != other.arr.size())
        throw std::invalid_argument("Size of expression doesn't match another.");

    if (exp_type != NO_SYMBOL || other.exp_type != NO_SYMBOL)
        throw std::invalid_argument("At least one of expressions has a symbol.");

    Expression result(*this);
    for (int i = 0; i < arr.size(); ++i) {
        result.arr[i] += other.arr[i];
    }

    result.b = b + other.b;

    return result;
}

Expression Expression::operator+(const double& scalar) const{
    Expression result(*this);
    result.b += scalar;

    return result;
}

Expression& Expression::operator+=(const LpVariable& var){
    if (arr.size() < var.id)
        throw std::invalid_argument("Size of expression doesn't match variables.");

    ++arr[var.id];

    return *this;
}

Expression& Expression::operator+=(const Expression& other){
    if (arr.size() != other.arr.size())
        throw std::invalid_argument("Size of expression doesn't match another.");

    if (exp_type != NO_SYMBOL || other.exp_type != NO_SYMBOL)
        throw std::invalid_argument("At least one of expressions has a symbol.");

    for (int i = 0; i < arr.size(); ++i) {
        arr[i] += other.arr[i];
    }
    b += other.b;

    return *this;
}

Expression& Expression::operator+=(const double& scalar){
    b += scalar;

    return *this;
}

Expression Expression::operator-(const LpVariable& var) const{
    if (arr.size() < var.id)
        throw std::invalid_argument("Size of expression doesn't match variables.");

    Expression result(*this);
    --result.arr[var.id];

    return result;
}

Expression Expression::operator-(const Expression& other) const{
    if (arr.size() != other.arr.size())
        throw std::invalid_argument("Size of expression doesn't match another.");

    if (exp_type != NO_SYMBOL || other.exp_type != NO_SYMBOL)
        throw std::invalid_argument("At least one of expressions has a symbol.");

    Expression result(*this);
    for (int i = 0; i < arr.size(); ++i) {
        result.arr[i] -= other.arr[i];
    }

    result.b -= other.b;

    return result;
}

Expression Expression::operator-(const double& scalar) const{
    Expression result(*this);
    result.b -= scalar;

    return result;
}

Expression& Expression::operator-=(const LpVariable& var){
    if (arr.size() < var.id)
        throw std::invalid_argument("Size of expression doesn't match variables.");

    --arr[var.id];

    return *this;
}

Expression& Expression::operator-=(const Expression& other){
    if (arr.size() != other.arr.size())
        throw std::invalid_argument("Size of expression doesn't match another.");

    if (exp_type != NO_SYMBOL || other.exp_type != NO_SYMBOL)
        throw std::invalid_argument("At least one of expressions has a symbol.");

    for (int i = 0; i < arr.size(); ++i) {
        arr[i] -= other.arr[i];
    }
    b -= other.b;

    return *this;
}

Expression& Expression::operator-=(const double& scalar){
    b -= scalar;

    return *this;
}

Expression Expression::operator*(const double& scalar) const{
    Expression result(*this);
    for (int i = 0; i < arr.size(); ++i)
        result.arr[i] *= scalar;
    result.b *= scalar;

    return result;
}

Expression& Expression::operator*=(const double& scalar){
    for (double& num: arr)
        num *= scalar;
    b *= scalar;

    return *this;
}

Expression Expression::operator/(const double& scalar) const{
    Expression result(*this);
    for (double& num: result.arr)
        num /= scalar;
    result.b /= scalar;

    return result;
}

Expression& Expression::operator/=(const double& scalar){
    for (double& num: arr)
        num /= scalar;
    b /= scalar;

    return *this;
}

Expression Expression::operator<=(const Expression& other) const{
    if (arr.size() != other.arr.size())
        throw std::invalid_argument("Size of expression doesn't match another.");

    if (exp_type != NO_SYMBOL || other.exp_type != NO_SYMBOL)
        throw std::invalid_argument("At least one of expressions has a symbol.");

    Expression result(*this);
    for (int i = 0; i < arr.size(); ++i)
        result.arr[i] -= other.arr[i];
    result.b = other.b - b;
    result.exp_type = LESS_EQUAL;

    return result;
}

Expression Expression::operator<=(const LpVariable& var) const{
    if (arr.size() < var.id)
        throw std::invalid_argument("Size of expression doesn't match variables.");

    Expression result(*this);

    --result.arr[var.id];
    result.b = -result.b;
    result.exp_type = LESS_EQUAL;

    return result;
}

Expression Expression::operator<=(const double& scalar) const{
    Expression result(*this);

    result.b = scalar - result.b;
    result.exp_type = LESS_EQUAL;

    return result;
}

Expression Expression::operator>=(const Expression& other) const{
    if (arr.size() != other.arr.size())
        throw std::invalid_argument("Size of expression doesn't match another.");

    if (exp_type != NO_SYMBOL || other.exp_type != NO_SYMBOL)
        throw std::invalid_argument("At least one of expressions has a symbol.");

    Expression result(arr.size(), false);
    for (int i = 0; i < arr.size(); ++i)
        result.arr[i] = other.arr[i] - arr[i];
    result.b -= other.b;
    result.exp_type = LESS_EQUAL;

    return result;
}

Expression Expression::operator>=(const LpVariable& var) const{
    if (arr.size() < var.id)
        throw std::invalid_argument("Size of expression doesn't match variables.");

    Expression result(*this);

    for (double& num: result.arr)
        num = -num;
    ++result.arr[var.id];
    result.exp_type = LESS_EQUAL;

    return result;
}

Expression Expression::operator>=(const double& scalar) const{
    Expression result(*this);

    for (double& num: result.arr)
        num = -num;
    result.b = result.b - scalar;
    result.exp_type = LESS_EQUAL;

    return result;
}

Expression Expression::operator==(const Expression& other) const{
    if (arr.size() != other.arr.size())
        throw std::invalid_argument("Size of expression doesn't match another.");

    if (exp_type != NO_SYMBOL || other.exp_type != NO_SYMBOL)
        throw std::invalid_argument("At least one of expressions has a symbol.");

    Expression result(*this);
    for (int i = 0; i < arr.size(); ++i)
        result.arr[i] -= other.arr[i];
    result.b = other.b - b;
    result.exp_type = EQUAL;

    return result;
}

Expression Expression::operator==(const LpVariable& var) const{
    if (arr.size() <= var.id)
        throw std::invalid_argument("Size of expression doesn't match variables.");

    Expression result(*this);

    --result.arr[var.id];
    result.b = -result.b;
    result.exp_type = EQUAL;

    return result;
}

Expression Expression::operator==(const double& scalar) const {
    Expression result(*this);

    result.b = scalar - result.b;
    result.exp_type = EQUAL;

    return result;
}

std::ostream& operator<<(std::ostream& os, const Expression& expr) {
    bool first_time = true;

    for (int i = 0; i < expr.arr.size(); ++i) {
        if (expr.arr[i] == 0)
            continue;

        if (first_time) {
            os << expr.arr[i] << " * x" << i;
            first_time = false;
        }
        else {
            if (expr.arr[i] > 0)
                os << " + " << expr.arr[i] << " * x" << i;
            else
                os << " - " << -expr.arr[i] << " * x" << i;
        }
    }

    if (expr.exp_type == Expression::NO_SYMBOL) {
        if (expr.b > 0)
            os << " + " << expr.b << " (NO_SYMBOL)";
        else if (expr.b < 0)
            os << " - " << -expr.b << " (NO_SYMBOL)";
    }
    else if (expr.exp_type == Expression::LESS_EQUAL)
        os << " <= " << expr.b;
    else
        os << " == " << expr.b;

    return os;
}

LpVariable::LpVariable(const double& low, const double& up) {
    if (low > up)
        throw std::invalid_argument("The upper bound must be greater than or equal to the lower bound.");

    lower_bound = low;
    upper_bound = up;
    value = 0;

    id = id_counter;
    ++id_counter;
}

Expression LpVariable::operator+(const LpVariable& other) const {
    Expression result(id_counter);

    result.arr[id] = 1;
    result.arr[other.id] = 1;

    return result;
}

Expression LpVariable::operator-(const LpVariable& other) const {
    Expression result(id_counter);

    result.arr[id] = 1;
    result.arr[other.id] = -1;

    return result;
}

Expression LpVariable::operator*(const double& scalar) const {
    Expression result(id_counter);

    result.arr[id] = scalar;

    return result;
}

Expression LpVariable::operator/(const double& scalar) const {
    Expression result(id_counter);

    result.arr[id] = 1 / scalar;

    return result;
}

void LpVariable::set_bound(const double &low, const double &up) {
    lower_bound = low > lower_bound ? low : lower_bound;
    upper_bound = up < upper_bound ? up : upper_bound;
}

void LpVariable::set_lower_bound(const double &low) {
    lower_bound = low > lower_bound ? low : lower_bound;
}

void LpVariable::set_upper_bound(const double &up) {
    upper_bound = up < upper_bound ? up : upper_bound;
}

std::ostream& operator<<(std::ostream& os, const LpVariable& var) {
    return os << 'x' << var.id;
}

Expression LpVariable::operator+(const Expression &expr) const {
    return expr + *this;
}

Expression LpVariable::operator-(const Expression &expr) const {
    if (expr.arr.size() < id)
        throw std::invalid_argument("Size of expression doesn't match variables.");

    Expression result(expr);
    for (int i = 0; i < id; ++i) {
        result.arr[i] = -result.arr[i];
    }
    for (int i = id + 1; i < result.arr.size(); ++i) {
        result.arr[i] = -result.arr[i];
    }
    result.arr[id] = 1 - result.arr[id];

    return result;
}

Expression operator*(const double &scalar, const LpVariable &var) {
    return var * scalar;
}
