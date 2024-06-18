#ifndef functions_LIBE_gm2_H
#define functions_LIBE_gm2_H
#include "non_linear_fit.hpp"
#include "functions_amu.hpp"

static constexpr double e_em = sqrt(4 * M_PI * alpha_em);

static constexpr double Mpip_exp = 139.57039;
static constexpr double Mpip_exp_err = 0.0001;
static constexpr double MKp_exp = 493.677;
static constexpr double MKp_exp_err = 0.016;
static constexpr double MK0_exp = 497.611;
static constexpr double MK0_exp_err = 0.013;

double lhs_function_LIBE_gm2_eg(int j, double**** in, int t, struct fit_type fit_info);

double lhs_M_correction(int j, double**** in, int t, struct fit_type fit_info);
double lhs_M_correction_fit(int j, double**** in, int t, struct fit_type fit_info);
double rhs_M_correction_fit(int n, int Nvar, double* x, int Npar, double* P);


double lhs_M_mefft_correction(int j, double**** in, int t, struct fit_type fit_info);

double lhs_Mt_correction(int j, double**** in, int t, struct fit_type fit_info);
double rhs_fit_mass_correction(int n, int Nvar, double* x, int Npar, double* P);

double** add_corr_sum_k_VKVK(int j, double**** in, int t, struct fit_type fit_info);
double** add_corr_correct_VKVK(int j, double**** in, int t, struct fit_type fit_info);
double** deriv_e(int j, double**** in, int t, struct fit_type fit_info);
double** deriv_e_exchange(int j, double**** in, int t, struct fit_type fit_info);
double** compute_TM_mass_insertion(int j, double**** in, int t, struct fit_type fit_info);
double** compute_TM_critical_mass_insertion(int j, double**** in, int t, struct fit_type fit_info);

double lhs_dm0_cr(int j, double**** in, int t, struct fit_type fit_info);
double lhs_dm0_cr_nabla(int j, double**** in, int t, struct fit_type fit_info);

struct fit_result  solve_QED_system(struct fit_type fit_info);

double rhs_fit_dmu_phys(int n, int Nvar, double* x, int Npar, double* P);


#endif
