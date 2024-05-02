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
    double r = (in[j][id_DC][(t + 1) % T][0] / in[j][id_C][(t + 1) % T][0]) - (in[j][id_DC][t][0] / in[j][id_C][t][0]);
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

    int ime1 = fit_info.corr_id[3];
    int ipe1 = fit_info.corr_id[5];

    int ime2 = fit_info.corr_id[1];
    int ipe2 = fit_info.corr_id[7];

    int imme = fit_info.corr_id[0];
    int ippe = fit_info.corr_id[8];

    int i00 = fit_info.corr_id[4];

    double e1 = fit_info.ave_P[0];
    double e2 = fit_info.ave_P[1];

    double e = fit_info.ave_P[2];

    double d1 = (in[j][ipe1][t][0] - 2 * in[j][i00][t][0] + in[j][ime1][t][0]) / (e * e);
    double d2 = (in[j][ipe2][t][0] - 2 * in[j][i00][t][0] + in[j][ime2][t][0]) / (e * e);
    double d3 = (in[j][ippe][t][0] - 2 * in[j][i00][t][0] + in[j][imme][t][0]) / (e * e);

    r[0][0] = (e1 * e2 * (d3 - d1 - d2) + e1 * e1 * d1 + e2 * e2 * d2) / 2.0;
    r[0][1] = 0;

    return r;
}
