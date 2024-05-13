#ifndef functions_LIBE_gm2_H
#define functions_LIBE_gm2_H
#include "non_linear_fit.hpp"

static constexpr double alpha_em = 1.0/137.035999174;
static constexpr double e_em = sqrt(4*M_PI*alpha_em);

double lhs_function_LIBE_gm2_eg(int j, double**** in, int t, struct fit_type fit_info);

double lhs_M_correction(int j, double**** in, int t, struct fit_type fit_info);
double lhs_M_mefft_correction(int j, double**** in, int t, struct fit_type fit_info);

double lhs_Mt_correction(int j, double**** in, int t, struct fit_type fit_info);
double rhs_fit_mass_correction(int n, int Nvar, double* x, int Npar, double* P);

double** deriv_e(int j, double**** in, int t, struct fit_type fit_info);
double** deriv_e_exchange(int j, double**** in, int t, struct fit_type fit_info);
double** compute_TM_mass_insertion(int j, double**** in, int t, struct fit_type fit_info);
double** compute_TM_critical_mass_insertion(int j, double**** in, int t, struct fit_type fit_info);

double lhs_dm0_cr(int j, double**** in, int t, struct fit_type fit_info);
double lhs_dm0_cr_nabla(int j, double**** in, int t, struct fit_type fit_info);
#endif
