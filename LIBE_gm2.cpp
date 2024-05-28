#define CONTROL
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// #include <time.h>
// #include <complex.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <map>

#include "global.hpp"
#include "read.hpp"
#include "resampling_new.hpp"
#include "tower.hpp"
#include "linear_fit.hpp"
#include "mutils.hpp"
#include "various_fits.hpp"

#include "correlators_analysis.hpp"
#include "fit_all.hpp"

#include "functions_LIBE_gm2.hpp"
struct kinematic kinematic_2pt;

generic_header read_head(FILE* stream) {
    generic_header header;
    return header;
}
void write_header_g2(FILE* jack_file, generic_header head) {
    int fi = 0;
    fi += fwrite(&head.T, sizeof(int), 1, jack_file);
    fi += fwrite(&head.L, sizeof(int), 1, jack_file);
    int nmu = head.mus.size();
    fi += fwrite(&nmu, sizeof(int), 1, jack_file);
    for (double mu : head.mus) {
        fi += fwrite(&mu, sizeof(double), 1, jack_file);
    }
}

char** argv_to_options(char** argv) {
    char** option;
    option = (char**)malloc(sizeof(char*) * 7);
    option[0] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[1] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[2] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[3] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[4] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[5] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[6] = (char*)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option[1], NAMESIZE, "read_plateaux"); // blind/see/read_plateaux
    mysprintf(option[2], NAMESIZE, "-p");            // -p
    mysprintf(option[3], NAMESIZE, argv[2]);         // path
    mysprintf(option[4], NAMESIZE, argv[6]);         // resampling
    mysprintf(option[5], NAMESIZE, "no");            // pdf
    mysprintf(option[6], NAMESIZE, argv[3]);         // infile
    return option;
}

void init_global_head(generic_header head) {
    file_head.l1 = head.L;
    file_head.l0 = head.T;
    file_head.l2 = head.L;
    file_head.l3 = head.L;
    file_head.nk = 2;
    file_head.musea = head.mus[1];
    file_head.k = (double*)malloc(sizeof(double) * file_head.nk * 2);
    file_head.k[0] = 0;
    file_head.k[1] = 0;
    file_head.k[2] = head.mus[1];
    file_head.k[3] = head.mus[1];

    file_head.nmoms = 1;
    file_head.mom = (double**)malloc(sizeof(double*) * file_head.nmoms);
    for (int i = 0; i < file_head.nmoms; i++) {
        file_head.mom[i] = (double*)malloc(sizeof(double) * 4);
        file_head.mom[i][0] = 0;
        file_head.mom[i][1] = 0;
        file_head.mom[i][2] = 0;
        file_head.mom[i][3] = 0;
    }
}

void read_twopt(FILE* stream, double*** to_write, generic_header head) {
    int fi = 0;
    int id;
    int i = fread(&id, sizeof(int), 1, stream);
    for (int k = 0; k < head.ncorr; k++) {
        for (int t = 0; t < head.T; t++) {
            fi += fread(to_write[k][t], sizeof(double), 2, stream);
        }
    }
}

int id_twpt(generic_header head, int im, int im1, int TMOS, int ig, int id_counter) {

    return id_counter + (head.oranges.size()) * (ig
        + head.gammas.size() * (TMOS + head.rs.size() * (im1 + 2 * (im))));
}

int main(int argc, char** argv) {
    error(argc != 7, 1, "nissa_mpcac ",
        "usage:./nissa_mpcac -p path file -bin $bin  jack/boot \n separate "
        "path and file please");

    char** option = argv_to_options(argv);

    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], option[6]);

    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");

    FILE* infile = open_file(namefile, "r");

    //////////////////////////////////// read and setup header
    generic_header head; //= read_head(infile);
    head.read_header_debug(infile);
    head.print_header();
    init_global_head(head);

    //////////////////////////////////// setup jackboot and binning
    int confs = head.Njack;
    int bin = atoi(argv[5]);
    int Neff = confs / bin;
    int Njack;
    if (strcmp(argv[6], "jack") == 0) {
        Njack = Neff + 1;
        myres = new resampling_jack(Neff);
    }
    else if (strcmp(argv[6], "boot") == 0) {
        Njack = (Neff * 2 + 1);
        myres = new resampling_boot(Neff * 2);
    }
    else {
        Njack = 0;
        error(1 == 1, 1, "main", "argv[7]= %s is not jack or boot", argv[7]);
    }
    // now Njack need to be the number of jacks
    head.Nconf = head.Njack;
    head.Njack = Njack;
    //////////////////////////////////// setup output files
    mysprintf(namefile, NAMESIZE, "%s/out/%s_output", option[3], option[6]);
    printf("writing output in :\n %s \n", namefile);
    FILE* outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/jackknife/%s_%s", option[3], option[4],
        option[6]);
    FILE* jack_file = open_file(namefile, "w+");
    // write_header_g2(jack_file, head);
    head.write_header(jack_file);

    //////////////////////////////////// confs
    double**** data = calloc_corr(confs, head.ncorr, head.T);

    printf("confs=%d\n", confs);
    printf("ncorr=%d\n", head.ncorr);
    printf("kappa=%g\n", head.kappa);
    for (int iconf = 0; iconf < confs; iconf++) {
        read_twopt(infile, data[iconf], head);
    }

    double**** data_bin = binning(confs, head.ncorr, head.T, data, bin);
    double**** conf_jack = myres->create(Neff, head.ncorr, head.T, data_bin);
    free_corr(Neff, head.ncorr, head.T, data_bin);
    free_corr(confs, head.ncorr, head.T, data);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // print all the effective masses correlators
    // set the option to not read for a plateaux
    mysprintf(namefile, NAMESIZE, "%s/out/%s_meff_correlators", option[3],
        option[6]);
    FILE* outfile_meff_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_raw_correlators", option[3],
        option[6]);
    FILE* outfile_raw_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_shifted_correlators", option[3],
        option[6]);
    FILE* outfile_shifted_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_log_meff_shifted", option[3],
        option[6]);
    FILE* outfile_log_meff_shifted = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_gamma", option[3], option[6]);
    FILE* out_gamma = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/%s_HLT_kernel", option[3],
        option[6]);
    FILE* outfile_HLT_kernel = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_HLT_AoverB", option[3],
        option[6]);
    FILE* outfile_HLT_AoverB = open_file(namefile, "w+");


    char save_option[NAMESIZE];
    sprintf(save_option, "%s", option[1]);
    sprintf(option[1], "blind");
    FILE* dev_null = open_file("/dev/null", "w");
    struct fit_type fit_info_silent;
    fit_info_silent.verbosity = -1;
    fit_info_silent.chi2_gap_jackboot = 1e+6;
    fit_info_silent.guess_per_jack = 0;

    for (int icorr = 0; icorr < head.ncorr; icorr++) {
        // log effective mass
        double* tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_meff_corr, icorr, "log", M_eff_log, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        // raw correlator
        file_head.l0 = head.T * 2;
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity_im,
            dev_null, fit_info_silent);
        free(tmp_meff_corr);
        file_head.l0 = head.T;
        // shifted correlator
        // tmp_meff_corr = plateau_correlator_function(
        //     option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        //     namefile_plateaux, outfile_shifted_corr, icorr, "shift_cor", shift_corr,
        //     dev_null, fit_info_silent);
        // free(tmp_meff_corr);
        // // log_meff shifted correlator
        // tmp_meff_corr = plateau_correlator_function(
        //     option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        //     namefile_plateaux, outfile_log_meff_shifted, icorr, "log_shift",
        //     M_eff_log_shift, dev_null, fit_info_silent);
        // free(tmp_meff_corr);
    }
    fit_info_silent.restore_default();
    sprintf(option[1], "%s", save_option); // restore option
    corr_counter = -1;

    ////////////////////////////////////////////////////////////
    // symmetrize
    ////////////////////////////////////////////////////////////
    int idTM = 0;
    int idOS = 1;
    int same_mass = 1;
    int diff_mass = 0;
    std::map<std::string, int> gamma_map{
        {"P5P5", 0},
        {"S0S0", 1},
        {"A1A1", 2},
        {"A2A2", 3},
        {"A3A3", 4},
        {"A0A0", 5},
        {"V1V1", 6},
        {"V2V2", 7},
        {"V3V3", 8},
        {"V0V0", 9},
        //
        {"P5P5", 10},
        {"S0P5", 11},
        {"A1P5", 12},
        {"A2P5", 13},
        {"A3P5", 14},
        {"A0P5", 15},
        {"V1P5", 16},
        {"V2P5", 17},
        {"V3P5", 18},
        {"V0P5", 19},
        //
        {"P5P5", 20},
        {"P5S0", 21},
        {"P5A1", 22},
        {"P5A2", 23},
        {"P5A3", 24},
        {"P5A0", 25},
        {"P5V1", 26},
        {"P5V2", 27},
        {"P5V3", 28},
        {"P5V0", 29}
    };
    std::map<std::string, int> counterterm_map{
        {"-e-e", 0},
        {"0-e", 1},
        {"e-e", 2},
        {"-e0", 3},
        {"00", 4},
        {"e0", 5},
        {"-ee", 6},
        {"0e", 7},
        {"ee", 8},
        {"00_sib", 9},
        {"01_sib", 10},
        {"10_sib", 11},
        {"02_sib", 12},
        {"20_sib", 13},
    };

    int id_PS = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["00"]);
    // symmetrise_jackboot(Njack, id_PS, head.T, conf_jack);

    int id_K = id_twpt(head, head.mus.size() - 1, diff_mass, idTM, gamma_map["P5P5"], counterterm_map["00"]);
    // symmetrise_jackboot(Njack, id_K, head.T, conf_jack);

    int id_K1 = id_twpt(head, head.mus.size() - 2, diff_mass, idTM, gamma_map["P5P5"], counterterm_map["00"]);
    // symmetrise_jackboot(Njack, id_K1, head.T, conf_jack);

    int id_K2 = id_twpt(head, head.mus.size() - 3, diff_mass, idTM, gamma_map["P5P5"], counterterm_map["00"]);
    // symmetrise_jackboot(Njack, id_K2, head.T, conf_jack);

    for (int i = 0;i < head.oranges.size();i++) {
        int j = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], i);
        symmetrise_jackboot(Njack, j, head.T, conf_jack);
        j = id_twpt(head, head.mus.size() - 1, diff_mass, idTM, gamma_map["P5P5"], i);
        symmetrise_jackboot(Njack, j, head.T, conf_jack);
        j = id_twpt(head, head.mus.size() - 2, diff_mass, idTM, gamma_map["P5P5"], i);
        symmetrise_jackboot(Njack, j, head.T, conf_jack);
        j = id_twpt(head, head.mus.size() - 3, diff_mass, idTM, gamma_map["P5P5"], i);
        symmetrise_jackboot(Njack, j, head.T, conf_jack);

        j = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["V0P5"], i);
        symmetrise_jackboot(Njack, j, head.T, conf_jack, -1);

        j = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5V0"], i);
        symmetrise_jackboot(Njack, j, head.T, conf_jack, -1);

    }

    // symmetrise_jackboot(Njack, 1, head.T, conf_jack, -1);

    ////////////////////////////////////////////////////////////
    // start fitting
    //////////////////////////////
    corr_counter = -1;

    double* M_PS = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, id_PS, "M_{PS}", M_eff_T, jack_file);
    check_correlatro_counter(0);

    double* M_K = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, id_K, "M_{K}", M_eff_T, jack_file);
    check_correlatro_counter(1);

    double* M_K1 = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, id_K1, "M_{K1}", M_eff_T, jack_file);
    check_correlatro_counter(2);

    double* M_K2 = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, id_K2, "M_{K2}", M_eff_T, jack_file);
    check_correlatro_counter(3);

    //////////////////////////////  Pi + //////////////////////////////////////////////////////////////
    int ncorr_new = head.ncorr;
    int id_de_pi;
    double* e_Mpi, * mu_Mpi, * m0_Mpi;
    double* e_MKp, * muu_MKp, * mus_MKp, * m0u_MKp, * m0s_MKp;
    double* e_MK0, * mud_MK0, * mus_MK0, * m0d_MK0, * m0s_MK0;
    double* dm0;

    {   ////////////////////////////////////////////////////////////
        //  derivative respect to e^2
        ////////////////////////////////////////////////////////////
        struct fit_type fit_info;
        fit_info.N = 1;
        fit_info.Njack = Njack;

        fit_info.T = head.T;
        // fit_info.myen = { TDs, TJW ,mu, nu };
        fit_info.n_ext_P = 0;
        // int id_pi_e = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["de"]);
        // int id_pi_e0 = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["no"]);
        // int id_pi_me = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["-de"]);
        // int id_pi_e = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["de"]);
        // int id_pi_e0 = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["no"]);
        // int id_pi_me = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["-de"]);
        std::vector<int> id_pi(head.oranges.size());
        for (int i = 0; i < head.oranges.size(); i++) {
            id_pi[i] = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], i);
        }

        fit_info.corr_id = id_pi;
        fit_info.myen = { 0 }; // real or imag
        fit_info.ave_P = { 0.6666666666666666, -0.3333333333333333 , head.oranges[2] };// u , bar d
        add_correlators(option, ncorr_new, conf_jack, deriv_e, fit_info);
        id_de_pi = ncorr_new - 1;
        struct fit_type tmp_info;
        tmp_info.codeplateaux = true;
        tmp_info.tmin = 2;
        tmp_info.tmax = 2;
        double* tmp = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_de_pi, "Delta_e_C_{PS}", identity, jack_file, tmp_info);
        free(tmp);
        check_correlatro_counter(4);
        fit_info.restore_default();
    }
    {
        struct fit_type fit_info;
        struct fit_result fit_out;

        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 1;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
        fit_info.ext_P[0] = M_PS;

        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.T = head.T;
        fit_info.corr_id = { id_de_pi , id_PS };
        fit_info.myen = { 1, 0 }; // sign , reim

        struct fit_result fit_me_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_e_M_{PS}", fit_info,
            jack_file);
        e_Mpi = myres->create_copy(fit_me_P5P5.P[0]);
        check_correlatro_counter(5);
        // free_fit_result(fit_info, fit_out);
        fit_info.restore_default();
        fit_me_P5P5.clear();
    }

    {   ////////////////////////////////////////////////////////////
        //  mass correction
        ////////////////////////////////////////////////////////////
        struct fit_type fit_info;


        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 1;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
        fit_info.ext_P[0] = M_PS;

        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.T = head.T;

        //////////////////////////////////////////  up
        fit_info.myen = { -1, 1 }; // sign , reim
        /// mu1, mu2 (the one daggered), insertion 1 and 3 are on mu2, insertion 2 and 4 are on mu1
        int id_dmu_u_pi = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], head.bananas[0] + 4);
        fit_info.corr_id = { id_dmu_u_pi , id_PS };

        struct fit_result fit_dmu_u_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_mu_u_M_{PS}", fit_info,
            jack_file);
        mu_Mpi = myres->create_copy(fit_dmu_u_P5P5.P[0]);
        // myres->mult(mu_Mpi, mu_Mpi, 2);
        check_correlatro_counter(6);

        // fit_dmu_P5P5.clear();

        //////////////////////////////////////////  down
        fit_info.myen = { +1,  1 };
        int id_dmu_d_pi = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], head.bananas[0] + 3);
        fit_info.corr_id = { id_dmu_d_pi , id_PS };

        struct fit_result fit_dmu_d_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_mu_d_M_{PS}", fit_info,
            jack_file);
        check_correlatro_counter(7);
        fit_info.restore_default();

    }


    {   ////////////////////////////////////////////////////////////
        //  critical mass correction
        ////////////////////////////////////////////////////////////
        struct fit_type fit_info;


        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 1;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
        fit_info.ext_P[0] = M_PS;

        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.T = head.T;

        //////////////////////////////////////////  up
        fit_info.myen = { -1, 0 };
        int id_dmu_u_pi = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], head.bananas[0] + 2);
        fit_info.corr_id = { id_dmu_u_pi , id_PS };

        struct fit_result fit_dmu_u_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_m0_u_M_{PS}", fit_info,
            jack_file);
        m0_Mpi = myres->create_copy(fit_dmu_u_P5P5.P[0]);
        myres->mult(m0_Mpi, m0_Mpi, 2);

        check_correlatro_counter(8);

        fit_dmu_u_P5P5.clear();

        //////////////////////////////////////////  down
        fit_info.myen = { -1, 0 };
        int id_dmu_d_pi = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], head.bananas[0] + 1);
        fit_info.corr_id = { id_dmu_d_pi , id_PS };

        struct fit_result fit_dmu_d_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_m0_d_M_{PS}", fit_info,
            jack_file);
        check_correlatro_counter(9);
        fit_info.restore_default();
        fit_dmu_d_P5P5.clear();


    }
    //////////////////////////////////////////////////////////////
    // exchange
    //////////////////////////////////////////////////////////////
    {   ////////////////////////////////////////////////////////////
        //  derivative respect to e^2
        ////////////////////////////////////////////////////////////
        struct fit_type fit_info;
        fit_info.N = 1;
        fit_info.Njack = Njack;

        fit_info.T = head.T;
        // fit_info.myen = { TDs, TJW ,mu, nu };
        fit_info.n_ext_P = 0;

        std::vector<int> id_pi(head.oranges.size());
        for (int i = 0; i < head.oranges.size(); i++) {
            id_pi[i] = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], i);
        }

        fit_info.corr_id = id_pi;
        fit_info.myen = { 0 }; // real or imag
        fit_info.ave_P = { 0.6666666666666666, -0.3333333333333333 , head.oranges[2] };// u , bar d
        add_correlators(option, ncorr_new, conf_jack, deriv_e_exchange, fit_info);
        int id_de_pi_exc = ncorr_new - 1;

        fit_info.restore_default();
        struct fit_result fit_out;

        fit_info.Nvar = 1;
        fit_info.Npar = 2;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 2;
        // fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
        fit_info.malloc_ext_P();
        for (int j = 0;j < fit_info.Njack;j++) {
            fit_info.ext_P[0][j] = M_PS[j];
            fit_info.ext_P[1][j] = head.T;
        }

        fit_info.function = rhs_fit_mass_correction;
        fit_info.linear_fit = true;
        fit_info.T = head.T;
        fit_info.corr_id = { id_de_pi_exc , id_PS };
        fit_info.myen = { 1, 0 }; // sign , reim
        // fit ratio of corr
        struct fit_result fit_me_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_Mt_correction, "Delta_e_exc_fit_M_{PS}", fit_info,
            jack_file);
        check_correlatro_counter(10);
        // fit_info.restore_default();
        fit_me_P5P5.clear();


        fit_info.Npar = 1;
        fit_info.function = constant_fit;
        fit_me_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_e_exc_M_{PS}", fit_info,
            jack_file);
        check_correlatro_counter(11);
        fit_me_P5P5.clear();

        fit_info.function = constant_fit;
        fit_me_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_mefft_correction, "Delta_e_exc_mefft_M_{PS}", fit_info,
            jack_file);
        check_correlatro_counter(12);
        fit_me_P5P5.clear();

        fit_info.restore_default();


    }

    ////////////////////////////////////////////////////////////////////////////////
    ///// critical mass mpcac
    ////////////////////////////////////////////////////////////////////////////////
    {
        int id_VP = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5V0"], counterterm_map["00"]);
        printf("id VP = %d\n", id_VP);
        struct fit_type tmp_info;
        tmp_info.codeplateaux = true;
        tmp_info.tmin = 2;
        tmp_info.tmax = 2;

        double* VP = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_VP, "VP", identity, jack_file, tmp_info);
        check_correlatro_counter(13);


        double* VP_im = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_VP, "VP_im", identity_im, jack_file, tmp_info);
        check_correlatro_counter(14);

        ////////////////////////////////////////////////////////////
        //  derivative respect to e^2
        ////////////////////////////////////////////////////////////
        struct fit_type fit_info;
        fit_info.N = 1;
        fit_info.Njack = Njack;

        fit_info.T = head.T;
        fit_info.n_ext_P = 0;
        std::vector<int> id_pi(head.oranges.size());
        for (int i = 0; i < head.oranges.size(); i++) {
            id_pi[i] = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5V0"], i);
        }

        fit_info.corr_id = id_pi;
        // fit_info.myen = { 1 }; // real or imag
        fit_info.ave_P = { 0.6666666666666666, -0.3333333333333333 , head.oranges[2] };// u , bar d
        add_correlators(option, ncorr_new, conf_jack, deriv_e, fit_info);
        int id_de_VP = ncorr_new - 1;

        double* tmp = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_de_VP, "Delta_e_VP", identity_im, jack_file, tmp_info);
        free(tmp);
        check_correlatro_counter(15);
        fit_info.restore_default();

        ////////////////////////////////////////////////////////////
        //  critical mass correction
        ////////////////////////////////////////////////////////////
        int id_VP_m0u = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5V0"], head.bananas[0] + 2);
        tmp = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_VP_m0u, "Delta_m0u_VP", identity_im, jack_file, tmp_info);
        check_correlatro_counter(16);
        free(tmp);

        int id_VP_m0d = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5V0"], head.bananas[0] + 1);
        tmp = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_VP_m0d, "Delta_m0d_VP", identity_im, jack_file, tmp_info);
        check_correlatro_counter(17);
        free(tmp);

        ////////////////////////////////////////////////////////////
        //  mass correction
        ////////////////////////////////////////////////////////////
        int id_VP_mu_u = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5V0"], head.bananas[0] + 4);
        tmp = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_VP_mu_u, "Delta_muu_VP", identity, jack_file, tmp_info);
        check_correlatro_counter(18);
        free(tmp);

        int id_VP_mu_d = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5V0"], head.bananas[0] + 3);
        tmp = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_VP_mu_d, "Delta_mud_VP", identity, jack_file, tmp_info);
        check_correlatro_counter(19);
        free(tmp);


        /////// find dm_cr
        fit_info.corr_id = { id_de_VP,  id_VP_m0u, id_VP_m0d };
        fit_info.myen = { 1, 1, 1 }; // re or im
        fit_info.ave_P = { e_em };
        fit_info.linear_fit = true;
        fit_info.T = head.T;
        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.function = constant_fit;
        fit_result dm0_cr = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_dm0_cr, "dm0_cr", fit_info,
            jack_file);
        check_correlatro_counter(20);
        dm0_cr.clear();

        dm0_cr = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_dm0_cr_nabla, "dm0_cr_nabla", fit_info,
            jack_file);
        check_correlatro_counter(21);
        dm0 = myres->create_copy(dm0_cr.P[0]);
    }


    /////////////////////////////  K + //////////////////////////////////////////////////////////////
    int id_de_Kp;

    {   ////////////////////////////////////////////////////////////
        //  derivative respect to e^2
        ////////////////////////////////////////////////////////////
        struct fit_type fit_info;
        fit_info.N = 1;
        fit_info.Njack = Njack;

        fit_info.T = head.T;
        // fit_info.myen = { TDs, TJW ,mu, nu };
        fit_info.n_ext_P = 0;
        // int id_pi_e = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["de"]);
        // int id_pi_e0 = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["no"]);
        // int id_pi_me = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["-de"]);
        // int id_pi_e = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["de"]);
        // int id_pi_e0 = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["no"]);
        // int id_pi_me = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["-de"]);
        std::vector<int> id_pi(head.oranges.size());
        for (int i = 0; i < head.oranges.size(); i++) {
            id_pi[i] = id_twpt(head, head.mus.size() - 1, diff_mass, idTM, gamma_map["P5P5"], i);
        }

        fit_info.corr_id = id_pi;
        fit_info.myen = { 0 }; // real or imag
        fit_info.ave_P = { 0.6666666666666666, -0.3333333333333333 , head.oranges[2] }; // u , bar s
        add_correlators(option, ncorr_new, conf_jack, deriv_e, fit_info);
        id_de_Kp = ncorr_new - 1;
        struct fit_type tmp_info;
        tmp_info.codeplateaux = true;
        tmp_info.tmin = 2;
        tmp_info.tmax = 2;
        double* tmp = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_de_Kp, "Delta_e_C_{Kp}", identity, jack_file, tmp_info);
        free(tmp);
        check_correlatro_counter(22);
        fit_info.restore_default();
    }
    {
        struct fit_type fit_info;
        struct fit_result fit_out;

        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 1;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
        fit_info.ext_P[0] = M_K;

        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.T = head.T;
        fit_info.corr_id = { id_de_Kp , id_K };
        fit_info.myen = { 1, 0 }; // sign , reim

        struct fit_result fit_me_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_e_M_{Kp}", fit_info,
            jack_file);
        check_correlatro_counter(23);
        e_MKp = myres->create_copy(fit_me_P5P5.P[0]);
        fit_me_P5P5.clear();
        fit_info.restore_default();
    }

    {   ////////////////////////////////////////////////////////////
        //  mass correction
        ////////////////////////////////////////////////////////////
        struct fit_type fit_info;


        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 1;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
        fit_info.ext_P[0] = M_K;

        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.T = head.T;

        //////////////////////////////////////////  up
        fit_info.myen = { -1, 1 }; // sign , reim
        int id_dmu_u_pi = id_twpt(head, head.mus.size() - 1, diff_mass, idTM, gamma_map["P5P5"], head.bananas[0] + 4);
        fit_info.corr_id = { id_dmu_u_pi , id_K };

        struct fit_result fit_dmu_u_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_mu_u_M_{Kp}", fit_info,
            jack_file);
        muu_MKp = myres->create_copy(fit_dmu_u_P5P5.P[0]);
        fit_dmu_u_P5P5.clear();
        check_correlatro_counter(24);

        //////////////////////////////////////////  down
        fit_info.myen = { +1,  1 };
        int id_dmu_d_pi = id_twpt(head, head.mus.size() - 1, diff_mass, idTM, gamma_map["P5P5"], head.bananas[0] + 3);
        fit_info.corr_id = { id_dmu_d_pi , id_K };

        struct fit_result fit_dmu_d_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_mu_d_M_{Kp}", fit_info,
            jack_file);
        check_correlatro_counter(25);
        mus_MKp = myres->create_copy(fit_dmu_d_P5P5.P[0]);
        fit_dmu_d_P5P5.clear();
        fit_info.restore_default();

    }


    {   ////////////////////////////////////////////////////////////
        //  critical mass correction
        ////////////////////////////////////////////////////////////
        struct fit_type fit_info;


        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 1;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
        fit_info.ext_P[0] = M_K;

        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.T = head.T;

        //////////////////////////////////////////  up
        fit_info.myen = { -1, 0 };
        int id_dmu_u_pi = id_twpt(head, head.mus.size() - 1, diff_mass, idTM, gamma_map["P5P5"], head.bananas[0] + 2);
        fit_info.corr_id = { id_dmu_u_pi , id_K };

        struct fit_result fit_dmu_u_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_m0_u_M_{Kp}", fit_info,
            jack_file);
        check_correlatro_counter(26);
        m0u_MKp = myres->create_copy(fit_dmu_u_P5P5.P[0]);
        fit_dmu_u_P5P5.clear();
        // fit_dmu_P5P5.clear();

        //////////////////////////////////////////  down
        fit_info.myen = { -1, 0 };
        int id_dmu_d_pi = id_twpt(head, head.mus.size() - 1, diff_mass, idTM, gamma_map["P5P5"], head.bananas[0] + 1);
        fit_info.corr_id = { id_dmu_d_pi , id_K };

        struct fit_result fit_dmu_d_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_m0_d_M_{Kp}", fit_info,
            jack_file);
        check_correlatro_counter(27);
        m0s_MKp = myres->create_copy(fit_dmu_d_P5P5.P[0]);
        fit_dmu_d_P5P5.clear();
        fit_info.restore_default();

    }

    /////////////////////////////  K0 //////////////////////////////////////////////////////////////
    int id_de_K0;

    {   ////////////////////////////////////////////////////////////
        //  derivative respect to e^2
        ////////////////////////////////////////////////////////////
        struct fit_type fit_info;
        fit_info.N = 1;
        fit_info.Njack = Njack;

        fit_info.T = head.T;
        // fit_info.myen = { TDs, TJW ,mu, nu };
        fit_info.n_ext_P = 0;
        // int id_pi_e = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["de"]);
        // int id_pi_e0 = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["no"]);
        // int id_pi_me = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["-de"]);
        // int id_pi_e = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["de"]);
        // int id_pi_e0 = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["no"]);
        // int id_pi_me = id_twpt(head, head.mus.size() - 1, same_mass, idTM, gamma_map["P5P5"], counterterm_map["-de"]);
        std::vector<int> id_pi(head.oranges.size());
        for (int i = 0; i < head.oranges.size(); i++) {
            id_pi[i] = id_twpt(head, head.mus.size() - 1, diff_mass, idTM, gamma_map["P5P5"], i);
        }

        fit_info.corr_id = id_pi;
        fit_info.myen = { 0 }; // real or imag
        fit_info.ave_P = { -0.3333333333333333, -0.3333333333333333 , head.oranges[2] };// d , bar s
        add_correlators(option, ncorr_new, conf_jack, deriv_e, fit_info);
        id_de_K0 = ncorr_new - 1;
        struct fit_type tmp_info;
        tmp_info.codeplateaux = true;
        tmp_info.tmin = 2;
        tmp_info.tmax = 2;
        double* tmp = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_de_K0, "Delta_e_C_{K0}", identity, jack_file, tmp_info);
        free(tmp);
        check_correlatro_counter(28);
        fit_info.restore_default();
    }
    {
        struct fit_type fit_info;
        struct fit_result fit_out;

        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 1;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
        fit_info.ext_P[0] = M_K;

        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.T = head.T;
        fit_info.corr_id = { id_de_K0 , id_K };
        fit_info.myen = { 1, 0 }; // sign , reim

        struct fit_result fit_me_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_e_M_{K0}", fit_info,
            jack_file);
        check_correlatro_counter(29);
        e_MK0 = myres->create_copy(fit_me_P5P5.P[0]);
        fit_me_P5P5.clear();

        // free_fit_result(fit_info, fit_out);
        fit_info.restore_default();
    }

    {   ////////////////////////////////////////////////////////////
        //  mass correction
        ////////////////////////////////////////////////////////////
        struct fit_type fit_info;


        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 1;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
        fit_info.ext_P[0] = M_K;

        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.T = head.T;

        //////////////////////////////////////////  up
        fit_info.myen = { -1, 1 }; // sign , reim
        int id_dmu_u_pi = id_twpt(head, head.mus.size() - 1, diff_mass, idTM, gamma_map["P5P5"], head.bananas[0] + 4);
        fit_info.corr_id = { id_dmu_u_pi , id_K };

        struct fit_result fit_dmu_u_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_mu_u_M_{K0}", fit_info,
            jack_file);
        mud_MK0 = myres->create_copy(fit_dmu_u_P5P5.P[0]);
        fit_dmu_u_P5P5.clear();

        check_correlatro_counter(30);

        // fit_dmu_P5P5.clear();

        //////////////////////////////////////////  down
        fit_info.myen = { +1,  1 };
        int id_dmu_d_pi = id_twpt(head, head.mus.size() - 1, diff_mass, idTM, gamma_map["P5P5"], head.bananas[0] + 3);
        fit_info.corr_id = { id_dmu_d_pi , id_K };

        struct fit_result fit_dmu_d_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_mu_d_M_{K0}", fit_info,
            jack_file);
        check_correlatro_counter(31);
        mus_MK0 = myres->create_copy(fit_dmu_d_P5P5.P[0]);
        fit_dmu_d_P5P5.clear();
        fit_info.restore_default();

    }


    {   ////////////////////////////////////////////////////////////
        //  critical mass correction
        ////////////////////////////////////////////////////////////
        struct fit_type fit_info;


        fit_info.Nvar = 1;
        fit_info.Npar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 1;
        fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);
        fit_info.ext_P[0] = M_K;

        fit_info.function = constant_fit;
        fit_info.linear_fit = true;
        fit_info.T = head.T;

        //////////////////////////////////////////  up
        fit_info.myen = { -1, 0 };
        int id_dmu_u_pi = id_twpt(head, head.mus.size() - 1, diff_mass, idTM, gamma_map["P5P5"], head.bananas[0] + 2);
        fit_info.corr_id = { id_dmu_u_pi , id_K };

        struct fit_result fit_dmu_u_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_m0_u_M_{K0}", fit_info,
            jack_file);
        m0d_MK0 = myres->create_copy(fit_dmu_u_P5P5.P[0]);
        fit_dmu_u_P5P5.clear();
        check_correlatro_counter(32);

        // fit_dmu_P5P5.clear();

        //////////////////////////////////////////  down
        fit_info.myen = { -1, 0 };
        int id_dmu_d_pi = id_twpt(head, head.mus.size() - 1, diff_mass, idTM, gamma_map["P5P5"], head.bananas[0] + 1);
        fit_info.corr_id = { id_dmu_d_pi , id_K };

        struct fit_result fit_dmu_d_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_m0_d_M_{K0}", fit_info,
            jack_file);
        m0s_MK0 = myres->create_copy(fit_dmu_d_P5P5.P[0]);
        fit_dmu_d_P5P5.clear();
        check_correlatro_counter(33);
        fit_info.restore_default();

    }

    //////////////////////////////////////////////////////////////
    // system
    // M * d = Mphys-Miso
    //////////////////////////////////////////////////////////////

    // create the jackknife
    double* Mpip_exp_jack = myres->create_fake(Mpip_exp, Mpip_exp_err, 1000);
    double* MKp_exp_jack = myres->create_fake(MKp_exp, MKp_exp_err, 1001);
    double* MK0_exp_jack = myres->create_fake(MK0_exp, MK0_exp_err, 1002);
    double mean, err;
    int seed;
    line_read_param(option, "a", mean, err, seed, namefile_plateaux);
    double* a_fm = myres->create_fake(mean, err, seed);
    double* a_MeV = myres->create_copy(a_fm);
    myres->div(a_MeV, a_fm, hbarc);
    {
        struct fit_type fit_info;
        int N = 3;
        fit_info.N = N;
        fit_info.Njack = Njack;
        double*** matrix = malloc_3<double>(Njack, N, N);
        double** b = malloc_2<double>(Njack, N);
        double** counter = malloc_2<double>(N, Njack);
        double** counter_Mev = malloc_2<double>(N, Njack);
        for (size_t j = 0; j < fit_info.Njack; j++) {
            matrix[j][0][0] = mu_Mpi[j];
            matrix[j][0][1] = mu_Mpi[j];
            matrix[j][0][2] = 0;

            matrix[j][1][0] = muu_MKp[j];
            matrix[j][1][1] = 0;
            matrix[j][1][2] = mus_MKp[j];

            matrix[j][2][0] = 0;
            matrix[j][2][1] = mud_MK0[j];
            matrix[j][2][2] = mus_MK0[j];

            b[j][0] = a_MeV[j] * Mpip_exp_jack[j] - M_PS[j] - e_em * e_em * e_Mpi[j] - 2 * dm0[j] * m0_Mpi[j];
            b[j][1] = a_MeV[j] * MKp_exp_jack[j] - M_K[j] - e_em * e_em * e_MKp[j] - (dm0[j] * 2) * (4.0 / 5.0) * m0u_MKp[j] - (dm0[j] * 2) * (1.0 / 5.0) * m0s_MKp[j];
            b[j][2] = a_MeV[j] * MK0_exp_jack[j] - M_K[j] - e_em * e_em * e_MK0[j] - (dm0[j] * 2) * (1.0 / 5.0) * m0d_MK0[j] - (dm0[j] * 2) * (1.0 / 5.0) * m0s_MK0[j];

            double* res = LU_decomposition_solver(N, matrix[j], b[j]);
            counter[0][j] = res[0];
            counter[1][j] = res[1];
            counter[2][j] = res[2];
            free(res);
            counter_Mev[0][j] = counter[0][j] / a_MeV[j];
            counter_Mev[1][j] = counter[1][j] / a_MeV[j];
            counter_Mev[2][j] = counter[2][j] / a_MeV[j];

        }
        printf("system:\n");
        double* err = (double*)malloc(sizeof(double) * Njack);
        for (int k = 0;k < N;k++) {
            for (int i = 0;i < N;i++) {
                for (size_t j = 0; j < fit_info.Njack; j++)
                    err[j] = matrix[j][k][i];
                printf("(%g  \\pm  %g )  ", matrix[Njack - 1][k][i], myres->comp_error(err));

            }
            printf("\n");
        }
        printf("*x=\n");
        for (int i = 0;i < N;i++) {

            for (size_t j = 0; j < fit_info.Njack; j++)
                err[j] = b[j][i];
            printf("(%g  \\pm  %g )  ", b[Njack - 1][i], myres->comp_error(err));
        }
        printf("\n");
        printf("e corrections  e^2= %g\n",e_em*e_em);
        printf("de_Mpi = (%g  \\pm  %g )  \n", e_Mpi[Njack - 1], myres->comp_error(e_Mpi));
        printf("de_MKp = (%g  \\pm  %g )  \n", e_MKp[Njack - 1], myres->comp_error(e_MKp));
        printf("de_MK0 = (%g  \\pm  %g )  \n", e_MK0[Njack - 1], myres->comp_error(e_MK0));
        printf("m0 corrections\n");
        printf("dm0_u = (%g  \\pm  %g )  \n", (dm0[Njack - 1] * 2) * (4.0 / 5.0), (myres->comp_error(dm0) * 2) * (4.0 / 5.0));
        printf("dm0_d = (%g  \\pm  %g )  \n", (dm0[Njack - 1] * 2) * (1.0 / 5.0), (myres->comp_error(dm0) * 2) * (4.0 / 5.0));
        printf("dm0u_Mpi = (%g  \\pm  %g )  \n", m0_Mpi[Njack - 1], myres->comp_error(m0_Mpi));
        printf("dm0d_Mpi = (%g  \\pm  %g )  \n", m0_Mpi[Njack - 1], myres->comp_error(m0_Mpi));
        printf("dm0u_MKp = (%g  \\pm  %g )  \n", m0u_MKp[Njack - 1], myres->comp_error(m0u_MKp));
        printf("dm0s_MKp = (%g  \\pm  %g )  \n", m0s_MKp[Njack - 1], myres->comp_error(m0s_MKp));
        printf("dm0d_MK0 = (%g  \\pm  %g )  \n", m0d_MK0[Njack - 1], myres->comp_error(m0d_MK0));
        printf("dm0s_MK0 = (%g  \\pm  %g )  \n", m0s_MK0[Njack - 1], myres->comp_error(m0s_MK0));
        free(err);

        printf("a \\delta \\mu_u = %g  \\pm  %g \\\\\n", counter[0][Njack - 1], myres->comp_error(counter[0]));
        printf("a \\delta \\mu_d = %g  \\pm %g  \\\\\n", counter[1][Njack - 1], myres->comp_error(counter[1]));
        printf("a \\delta \\mu_s = %g  \\pm %g  \\\\\n", counter[2][Njack - 1], myres->comp_error(counter[2]));
        printf("\\delta \\mu_u[Mev] = %g \\pm  %g  \\\\\n", counter_Mev[0][Njack - 1], myres->comp_error(counter_Mev[0]));
        printf("\\delta \\mu_d[Mev] = %g \\pm  %g  \\\\\n", counter_Mev[1][Njack - 1], myres->comp_error(counter_Mev[1]));
        printf("\\delta \\mu_s[Mev] = %g \\pm  %g  \\\\\n", counter_Mev[2][Njack - 1], myres->comp_error(counter_Mev[2]));

    }
    // free stuff
    free(M_PS);
    free(M_K);
    free(M_K1);
    free(M_K2);
    for (int i = 0; i < 7; i++) free(option[i]);
    free(option);
}