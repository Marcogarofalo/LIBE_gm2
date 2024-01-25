#ifndef functions_LIBE_gm2_H
#define functions_LIBE_gm2_H
#include "non_linear_fit.hpp"

double lhs_function_LIBE_gm2_eg(int j, double**** in, int t, struct fit_type fit_info);

double lhs_M_correction(int j, double**** in, int t, struct fit_type fit_info);
double** deriv_e(int j, double**** in, int t, struct fit_type fit_info);

#endif
