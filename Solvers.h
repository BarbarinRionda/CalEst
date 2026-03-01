// Solvers.h

#ifndef SOLVERS_H
#define SOLVERS_H

class GaussJordan {
public:
    GaussJordan();
    void solve();
};

class FactCrout {
public:
    FactCrout();
    void factorize();
};

#endif // SOLVERS_H
