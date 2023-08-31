#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <math.h>
#include <stdbool.h>

typedef double DTYPE_t;
typedef long long int ITYPE_t;

// Inline functions for Maximum Likelihood and Bayesian Information Criteria
DTYPE_t mll_mean(DTYPE_t x, DTYPE_t x2, DTYPE_t x3, ITYPE_t n) {
    return (x2 - (x * x) * 1.0 / (n));
}

DTYPE_t mll_var(DTYPE_t x, DTYPE_t x2, DTYPE_t x3, ITYPE_t n) {
    return (n * (log(2 * M_PI) + log(fmax(x3, 0.000000001) / n) + 1));
}


// Function for non-parametric Maximum Likelihood
DTYPE_t mll_nonparametric_ed(DTYPE_t* sumstatout, ITYPE_t tstar, ITYPE_t checklist, ITYPE_t nquantiles, ITYPE_t n) {
    DTYPE_t Fkl, temp_cost, cost = 0;
    ITYPE_t nseg = tstar - checklist;
    for (ITYPE_t isum = 0; isum < nquantiles; ++isum) {
        Fkl = sumstatout[isum] / nseg;
        temp_cost = (tstar - checklist) * (Fkl * log(Fkl) + (1 - Fkl) * log(1 - Fkl));
        cost += isnan(temp_cost) ? 0 : temp_cost;
    }
    return -2 * log(2 * n - 1) * cost / nquantiles;
}

DTYPE_t mll_meanvar(DTYPE_t x, DTYPE_t x2, DTYPE_t x3, ITYPE_t n) {
    return (n * (log(2 * M_PI) + log(fmax((x2 - ((x * x) / n)) / n, 0.00000000001)) + 1));
}

DTYPE_t mll_meanvar_exp(DTYPE_t x, DTYPE_t x2, DTYPE_t x3, ITYPE_t n) {
    return 2 * n * (log(x) - log(n));
}

DTYPE_t mll_meanvar_poisson(DTYPE_t x, DTYPE_t x2, DTYPE_t x3, ITYPE_t n) {
    return x == 0 ? 0 : 2 * x * (log(n) - log(x));
}

DTYPE_t mbic_var(DTYPE_t x, DTYPE_t x2, DTYPE_t x3, ITYPE_t n) {
    return (n * (log(2 * M_PI) + log(fmax(x3, 0.00000000001) / n) + 1) + log(n));
}

DTYPE_t mbic_meanvar(DTYPE_t x, DTYPE_t x2, DTYPE_t x3, ITYPE_t n) {
    return (n * (log(2 * M_PI) + log(fmax((x2 - ((x * x) / n)) / n, 0.00000000001)) + 1) + log(n));
}

DTYPE_t mbic_mean(DTYPE_t x, DTYPE_t x2, DTYPE_t x3, ITYPE_t n) {
    return (x2 - (x * x) / n + log(n));
}

DTYPE_t mbic_meanvar_exp(DTYPE_t x, DTYPE_t x2, DTYPE_t x3, ITYPE_t n) {
    return (2 * n * (log(x) - log(n)) + log(n));
}

DTYPE_t mbic_meanvar_poisson(DTYPE_t x, DTYPE_t x2, DTYPE_t x3, ITYPE_t n) {
    return x == 0 ? 0 : (2 * x * (log(n) - log(x)) + log(n));
}

DTYPE_t mll_nonparametric_ed_mbic(DTYPE_t* sumstatout, ITYPE_t tstar, ITYPE_t checklist, ITYPE_t nquantiles, ITYPE_t n) {
    DTYPE_t Fkl, temp_cost, cost = 0;
    ITYPE_t nseg = tstar - checklist;
    for (ITYPE_t isum = 0; isum < nquantiles; ++isum) {
        Fkl = sumstatout[isum] / nseg;
        temp_cost = (tstar - checklist) * (Fkl * log(Fkl) + (1 - Fkl) * log(1 - Fkl));
        cost += isnan(temp_cost) ? 0 : temp_cost;
    }
    return -2 * log(2 * n - 1) * cost / nquantiles;
}
