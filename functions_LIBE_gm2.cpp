#define functions_LIBE_gm2_C
#include "functions_LIBE_gm2.hpp"
#include "tower.hpp"

double lhs_function_LIBE_gm2_eg(int j, double**** in, int t, struct fit_type fit_info) {
    double r = in[j][0][t][0];
    return r;
}

double lhs_M_correction(int j, double**** in, int t, struct fit_type fit_info) {
    int id_DC = fit_info.corr_id[0];
    int id_C = fit_info.corr_id[1];
    int T = fit_info.T;
    double M = fit_info.ext_P[0][j];

    double den = (T / 2 - (t + 1)) * tanh(M * (T / 2 - (t + 1)));
    den -= (T / 2 - (t)) * tanh(M * (T / 2 - (t)));
    double r = (in[j][id_DC][(t + 1) % T][0] - in[j][id_DC][t][0]) / in[j][id_C][t][0];
    // if(j==fit_info.Njack-1){
    // printf("DC(%d)=%g\n", t, in[j][id_DC][t][0]);
    // printf("C(%d)=%g   %d\n", t, in[j][id_C][t][0],id_C);
    // printf("partial_DC(%d)=%g\n", t, in[j][id_DC][(t + 1) % T][0] - in[j][id_DC][t][0]);
    // printf("den(%d)=%g  %g\n", t,  tanh(M * (T / 2 - (t + 1))),M);
    // }
    r /= den;
    return r;
}


double** deriv_e(int j, double**** in, int t, struct fit_type fit_info) {
    double** r = malloc_2<double>(1, 2);
    int ie = fit_info.corr_id[0];
    int i0 = fit_info.corr_id[1];
    int ime = fit_info.corr_id[2];


    r[0][0] = (in[j][ie][t][0] - 2 * in[j][i0][t][0] + in[j][ime][t][0]) / 2.0;
    r[0][1] = 0;

    return r;

}
