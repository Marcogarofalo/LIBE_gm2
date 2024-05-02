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
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_shifted_corr, icorr, "shift_cor", shift_corr,
            dev_null, fit_info_silent);
        free(tmp_meff_corr);
        // log_meff shifted correlator
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_log_meff_shifted, icorr, "log_shift",
            M_eff_log_shift, dev_null, fit_info_silent);
        free(tmp_meff_corr);
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
        {"A4A4", 5},
        {"V1V1", 6},
        {"V2V2", 7},
        {"V3V3", 8},
        {"V4V4", 9}
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
        {"10_sib", 10},
        {"01_sib", 11},
        {"20_sib", 12},
        {"02_sib", 13},
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

    ////////////////////////////////////////////////////////////////////////////////////////////
    int ncorr_new = head.ncorr;
    int id_de_pi;

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
        fit_info.ave_P = { 0.6666666666666666, -0.3333333333333333 , head.oranges[2] };
        add_correlators(option, ncorr_new, conf_jack, deriv_e, fit_info);
        id_de_pi = ncorr_new - 1;
        struct fit_type tmp_info;
        tmp_info.codeplateaux = true;
        tmp_info.tmin = 2;
        tmp_info.tmax = 2;
        double* tmp = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_de_pi, "Delta_C_{PS}", identity, jack_file, tmp_info);
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

        struct fit_result fit_me_P5P5 = fit_fun_to_fun_of_corr(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
            outfile, lhs_M_correction, "Delta_M_{PS}", fit_info,
            jack_file);
        check_correlatro_counter(5);
        // free_fit_result(fit_info, fit_out);
        fit_info.restore_default();
        fit_me_P5P5.clear();
    }

    {   ////////////////////////////////////////////////////////////
        //  mass correction
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
        fit_info.ave_P = { 0.6666666666666666, -0.3333333333333333 , head.oranges[2] };
        add_correlators(option, ncorr_new, conf_jack, deriv_e, fit_info);
        id_de_pi = ncorr_new - 1;
        struct fit_type tmp_info;
        tmp_info.codeplateaux = true;
        tmp_info.tmin = 2;
        tmp_info.tmax = 2;
        double* tmp = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, id_de_pi, "Delta_C_{PS}", identity, jack_file, tmp_info);
        free(tmp);
        check_correlatro_counter(4);
        fit_info.restore_default();
    }

    // free stuff
    free(M_PS);
    free(M_K);
    free(M_K1);
    free(M_K2);
    for (int i = 0; i < 7; i++) free(option[i]);
    free(option);
}