#ifndef SOLVERS_H
#define SOLVERS_H

#include <vector>

class GaussJordan {
public:
    static std::vector<double> solve(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
};

class FactCrout {
public:
    static std::vector<double> solve(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
};

#endif // SOLVERS_H