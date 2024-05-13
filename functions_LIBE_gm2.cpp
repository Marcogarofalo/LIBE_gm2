#define functions_LIBE_gm2_C
#include "functions_LIBE_gm2.hpp"
#include "tower.hpp"
#include "correlators_analysis.hpp"

double lhs_function_LIBE_gm2_eg(int j, double**** in, int t, struct fit_type fit_info) {
    double r = in[j][0][t][0];
    return r;
}

double lhs_M_correction(int j, double**** in, int t, struct fit_type fit_info) {
    int id_DC = fit_info.corr_id[0];
    int id_C = fit_info.corr_id[1];
    int T = fit_info.T;
    double M = fit_info.ext_P[0][j];
    int  sign = fit_info.myen[0];
    int  reim = fit_info.myen[1];

    double den = (T / 2 - (t + 1)) * tanh(M * (T / 2 - (t + 1)));
    den -= (T / 2 - (t)) * tanh(M * (T / 2 - (t)));
    double r = (in[j][id_DC][(t + 1) % T][reim] / in[j][id_C][(t + 1) % T][0]) - (in[j][id_DC][t][reim] / in[j][id_C][t][0]);
    // if(j==fit_info.Njack-1){
    // printf("DC(%d)=%g\n", t, in[j][id_DC][t][0]);
    // printf("C(%d)=%g   %d\n", t, in[j][id_C][t][0],id_C);
    // printf("partial_DC(%d)=%g\n", t, in[j][id_DC][(t + 1) % T][0] - in[j][id_DC][t][0]);
    // printf("den(%d)=%g  %g\n", t,  tanh(M * (T / 2 - (t + 1))),M);
    // }
    r /= den;
    r *= sign;
    return r;
}


double lhs_Mt_correction(int j, double**** in, int t, struct fit_type fit_info) {
    int id_DC = fit_info.corr_id[0];
    int id_C = fit_info.corr_id[1];
    int T = fit_info.T;
    double M = fit_info.ext_P[0][j];
    int  sign = fit_info.myen[0];
    int  reim = fit_info.myen[1];

    double r = in[j][id_DC][t][reim] / in[j][id_C][t][0];

    r *= sign;
    return r;
}

double lhs_M_mefft_correction(int j, double**** in, int t, struct fit_type fit_info) {
    int id_DC = fit_info.corr_id[0];
    int id_C = fit_info.corr_id[1];
    int T = fit_info.T;
    double ct = in[j][id_C][t][0];
    double ctp1 = in[j][id_C][(t + 1 % T)][0];
    double M = M_eff_T_ct_ctp1(t, T, ct, ctp1);
    int  sign = fit_info.myen[0];
    int  reim = fit_info.myen[1];

    double den = (T / 2 - (t + 1)) * tanh(M * (T / 2 - (t + 1)));
    den -= (T / 2 - (t)) * tanh(M * (T / 2 - (t)));
    double r = (in[j][id_DC][(t + 1) % T][reim] / in[j][id_C][(t + 1) % T][0]) - (in[j][id_DC][t][reim] / in[j][id_C][t][0]);
    r /= den;
    r *= sign;
    return r;
}
double lhs_Mt_mefft_correction(int j, double**** in, int t, struct fit_type fit_info) {
    int id_DC = fit_info.corr_id[0];
    int id_C = fit_info.corr_id[1];
    int T = fit_info.T;
    int  sign = fit_info.myen[0];
    int  reim = fit_info.myen[1];


    double r = in[j][id_DC][t][reim] / in[j][id_C][t][0];

    r *= sign;
    return r;
}


double rhs_fit_mass_correction(int n, int Nvar, double* x, int Npar, double* P) {
    double t = x[0];
    int T = x[2];
    double M = x[1];
    double factor = (T / 2 - (t)) * tanh(M * (T / 2 - (t)));
    double r = P[1] + P[0] * factor;

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
    for (int reim = 0;reim < 2;reim++) {
        double d1 = (in[j][ipe1][t][reim] - 2 * in[j][i00][t][reim] + in[j][ime1][t][reim]) / (e * e);
        double d2 = (in[j][ipe2][t][reim] - 2 * in[j][i00][t][reim] + in[j][ime2][t][reim]) / (e * e);
        double d3 = (in[j][ippe][t][reim] - 2 * in[j][i00][t][reim] + in[j][imme][t][reim]) / (e * e);

        r[0][reim] = (e1 * e2 * (d3 - d1 - d2) + e1 * e1 * d1 + e2 * e2 * d2) / 2.0;

    }

    return r;
}

double** deriv_e_exchange(int j, double**** in, int t, struct fit_type fit_info) {

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
    for (int reim = 0;reim < 2;reim++) {

        double d1 = (in[j][ipe1][t][reim] - 2 * in[j][i00][t][reim] + in[j][ime1][t][reim]) / (e * e);
        double d2 = (in[j][ipe2][t][reim] - 2 * in[j][i00][t][reim] + in[j][ime2][t][reim]) / (e * e);
        double d3 = (in[j][ippe][t][reim] - 2 * in[j][i00][t][reim] + in[j][imme][t][reim]) / (e * e);

        r[0][reim] = (e1 * e2 * (d3 - d1 - d2)) / 2.0;
    }

    return r;
}


double** compute_TM_mass_insertion(int j, double**** in, int t, struct fit_type fit_info) {
    double** r = malloc_2<double>(1, 2);

    int im1 = fit_info.corr_id[0];
    int im2 = fit_info.corr_id[1];

    r[0][0] = -in[j][im1][t][1] + in[j][im2][t][1];
    r[0][1] = 0;

    return r;
}

double** compute_TM_critical_mass_insertion(int j, double**** in, int t, struct fit_type fit_info) {
    double** r = malloc_2<double>(1, 2);

    int im1 = fit_info.corr_id[0];
    int im2 = fit_info.corr_id[1];

    r[0][0] = -in[j][im1][t][0] - in[j][im2][t][0];
    r[0][1] = 0;

    return r;
}


double lhs_dm0_cr(int j, double**** in, int t, struct fit_type fit_info) {
    int id_e = fit_info.corr_id[0];
    int id_m0u = fit_info.corr_id[1];
    int id_m0d = fit_info.corr_id[2];
    int T = fit_info.T;
    int reim_e = fit_info.myen[0];
    int reim_m0u = fit_info.myen[1];
    int reim_m0d = fit_info.myen[2];
    double e = fit_info.ave_P[0];


    double r = e * e * in[j][id_e][t][reim_e] / (-in[j][id_m0u][t][reim_m0u] - in[j][id_m0d][t][reim_m0d]);

    return r;
}

double lhs_dm0_cr_nabla(int j, double**** in, int t, struct fit_type fit_info) {
    int id_e = fit_info.corr_id[0];
    int id_m0u = fit_info.corr_id[1];
    int id_m0d = fit_info.corr_id[2];
    int T = fit_info.T;
    int reim_e = fit_info.myen[0];
    int reim_m0u = fit_info.myen[1];
    int reim_m0d = fit_info.myen[2];
    double e = fit_info.ave_P[0];

    double num = in[j][id_e][(t + 1) % T][reim_e] - in[j][id_e][t][reim_e];
    double den = in[j][id_m0u][(t + 1) % T][reim_m0u] - in[j][id_m0u][t][reim_m0u];
    den += in[j][id_m0d][(t + 1) % T][reim_m0d] - in[j][id_m0d][t][reim_m0d];
    double r = e * e * num / (-den);

    return r;
}
