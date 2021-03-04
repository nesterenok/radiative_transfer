//
#include <omp.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "constants.h"
#include "spectroscopy.h"
#include "coll_rates_ch3oh.h"
#include "coll_rates_nh3.h"
#include "coll_rates_oh.h"
#include "coll_rates_h2o.h"
#include "coll_rates_h2co.h"

#include "cloud_data.h"
#include "iteration_lvg.h"
#include "iteration_control.h"
#include "transition_data.h"
#include "transition_list.h"
#include "maser_luminosity.h"

#define SOURCE_NAME "radiative_transfer.cpp"
#define NB_JOINED_LAYERS 2
#define MAX_NB_ITER_EXT 15000
#define MAX_NB_ITER_ACC 150

// These parameters must be the same as in the C-type shock simulations (are necessary for dust model),
#define HE_TO_H_NB_RATIO 0.09
#define STANDARD_NB_CR_PHOTONS 3000.
#define MICROTURBULENT_SPEED 3.e+4  // in cm/s, must be the same as in the shock simulations,
#define BEAMING_FACTOR 0.001  // = dOmega/4pi, for a cylindrical maser spot with line of sight length = 10 spot size, dOmega = 0.01 ster

using namespace std;

const double a_meth_fraction = 0.5; // A-methanol fraction,
const double oh2o_fraction = 0.75;  // ortho-H2O fraction, Emprechtinger et al., ApJ 765, p. 61 (2013); 
const double onh3_fraction = 0.5;   // ortho-NH3 fraction,
const double oh2co_fraction = 0.75; // ortho-H2CO fraction,


// the relative population error, is used in calculations of populations and in finding of inverted transitions,
const double rel_population_error = 1.e-5;
//const double population_error = 1.e-8; // not used

// OH radiative transfer calculations - with line overlap, population densities are saved
void calc_oh_masers(string input_data_path, string sim_data_path, string output_path, bool save_line_stat, int verbosity);

// OH radiative transfer calculations - no line overlap,
void calc_oh_masers_test(string input_data_path, string sim_data_path, string output_path, int verbosity);

// CH3OH radiative transfer calculations,
void calc_ch3oh_a_masers(string input_data_path, string sim_data_path, string output_path, bool save_lev_pop, int verbosity);
void calc_ch3oh_e_masers(string input_data_path, string sim_data_path, string output_path, int verbosity);

// H2O radiative transfer calculations, population densities are saved
// the data on population densities obtained by shock modeling may be used (for test)
void calc_h2o_para_masers(string input_data_path, string sim_data_path, string output_path, bool save_lev_pop, bool is_shock_data_used, int verbosity);
void calc_h2o_ortho_masers(string input_data_path, string sim_data_path, string output_path, bool save_lev_pop, int verbosity);

// NH3 radiative transfer calculations,
void calc_nh3_para_masers(string input_data_path, string sim_data_path, string output_path, int verbosity);
void calc_nh3_ortho_masers(string input_data_path, string sim_data_path, string output_path, int verbosity);

// H2CO radiative transfer calculations,
void calc_h2co_para_masers(string input_data_path, string sim_data_path, string output_path, int verbosity);
void calc_h2co_ortho_masers(string input_data_path, string sim_data_path, string output_path, int verbosity);


// the calculation of the line overlap number as a function of doppler width,
void calc_nb_line_overlaps(string input_data_path, string output_path);

// CH3OH collisions extrapolation
void ch3oh_coll_extrapolations(string input_data_path, string output_path);

void calc_molecular_populations(cloud_data* cloud, iteration_scheme_lvg* it_scheme_lvg, 
    energy_diagram* mol_levels, einstein_coeff* mol_einst, collisional_transitions* mol_coll, double* mol_popul, int nb_lev, 
    bool acceleration, int verbosity);

// saving the level populations in the format of the c-type shock simulations
void save_mol_data(std::string fname, cloud_data* cloud, double *lev_pop, int nb_lev);
void save_mol_vibr_data(std::string fname, cloud_data* cloud, energy_diagram* mol_levels, double* lev_pop, int nb_lev);  // free format


int main()
{
#ifdef _OPENMP
    omp_set_num_threads(4);
#pragma omp parallel 
    {
#pragma omp master 
        {
            cout << "OpenMP is supported" << endl;
            cout << "Nb of threads: " << omp_get_num_threads() << endl;
        }
    }
#endif

    bool save_line_stat, save_lev_pop, is_shock_data_used;
    int verbosity(1);
    string path, suff1, suff2, sim_data_path, input_data_path, output_path;
    
    suff1 = "2e6";
    suff2 = "1-16";

    // the path to the directory with spectroscopic data,
#ifdef __linux__
    input_data_path = "/disk2/nester/input_data/";
    path = "/disk2/nester/sim_data_2020_07/";

    stringstream lin_out;
    lin_out.clear();
    lin_out.str("");
    lin_out << "output_" << suff1 << "_" << suff2 << ".txt";

    ofstream outerr(lin_out.str().c_str(), ios::app);
    streambuf* orig_cerr = cerr.rdbuf(outerr.rdbuf());

    ofstream out(lin_out.str().c_str(), ios::app);
    streambuf* orig_cout = cout.rdbuf(out.rdbuf());
#else
    input_data_path = "C:/Users/Александр/Documents/input_data/";
    path = "C:/Users/Александр/Documents/Данные и графики/paper Cosmic masers in C-type shocks/";
#endif

    // calculation of OH line overlap statistics;
//    calc_nb_line_overlaps(input_data_path, output_path = "");
//    ch3oh_coll_extrapolations(input_data_path, output_path = "");

// check - extrapolated or not data on CH3OH collisions
    sim_data_path = path + "output_data_2e5/shock_cr3-15_17-5/";    
//    output_path = path + "ch3oh_2e5_extcollcoef/shock_cr1-16_17-5/";
//    calc_ch3oh_a_masers(input_data_path, sim_data_path, output_path, save_lev_pop = true, verbosity);
//    calc_ch3oh_e_masers(input_data_path, sim_data_path, output_path, verbosity);

    output_path = path + "oh_2e5/shock_cr3-15_17-5/";
    calc_oh_masers(input_data_path, sim_data_path, output_path, save_line_stat = true, verbosity); // with line statistics
 /*
//    output_path = "C:/Users/Александр/Documents/Данные и графики/paper Cosmic masers in C-type shocks/oh_nolo_2e5/shock_cr1-16_15/";
//    calc_oh_masers_test(input_data_path, sim_data_path, output_path, verbosity);

    save_line_stat = false;
    save_lev_pop = false; 
    is_shock_data_used = false;
    verbosity = false;
#pragma omp parallel shared(input_data_path, path, suff1, suff2, save_line_stat, save_lev_pop, is_shock_data_used, verbosity), private(output_path, sim_data_path)
    {
        int i, j, k;
        double time_in, time_tot, v;
        string suff, suff3;
        stringstream sstr;

#pragma omp for schedule (dynamic, 1)
        for (k = 0; k <= 40; k += 1) 
        {
            time_in = omp_get_wtime();
            sstr.clear();
            sstr.str("");

            v = 5. + 2.5 * k;
            i = (int)(v);
            j = (int)(v * 10);

            if (v < 9.999)
                sstr << "0";
            sstr << i;

            if (j > 10 * i)
                sstr << "-" << (j - i * 10);
            suff3 = sstr.str();

            suff = suff1 + "/shock_cr";
            suff += suff2 + "_";
            suff += suff3 + "/";
          
            sim_data_path = path + "output_data_";
            sim_data_path += suff;
            cout << "Simulation data: " << sim_data_path << endl;

            // line overlap is taken into account,
            // check - extended or usual data on OH-H2 collisions
            output_path = path + "oh_";
            output_path += suff;
            calc_oh_masers(input_data_path, sim_data_path, output_path, save_line_stat, verbosity);  // no line statistics

            // check - extrapolated or not data on CH3OH collisions
            output_path = path + "ch3oh_";
            output_path += suff;
            calc_ch3oh_a_masers(input_data_path, sim_data_path, output_path, save_lev_pop, verbosity);
            calc_ch3oh_e_masers(input_data_path, sim_data_path, output_path, verbosity);
                        
            output_path = path + "h2o_";
            output_path += suff;
            //calc_h2o_para_masers(input_data_path, sim_data_path, output_path, save_lev_pop, is_shock_data_used, verbosity);
            //calc_h2o_ortho_masers(input_data_path, sim_data_path, output_path, save_lev_pop, verbosity);

//            output_path = path + "nh3_";
//            output_path += suff;
//            calc_nh3_para_masers(input_data_path, sim_data_path, output_path, verbosity);
//            calc_nh3_ortho_masers(input_data_path, sim_data_path, output_path, verbosity);
            
//            output_path = path + "h2co_";
//            output_path += suff;
//            calc_h2co_para_masers(input_data_path, sim_data_path, output_path, verbosity);
//            calc_h2co_ortho_masers(input_data_path, sim_data_path, output_path, verbosity);
           
            time_tot = omp_get_wtime() - time_in;
            cout << "time in s for " << omp_get_thread_num() << ": " << time_tot << endl;
        }
    }*/
}

void calc_molecular_populations(cloud_data* cloud, iteration_scheme_lvg* it_scheme_lvg, 
    energy_diagram *mol_levels, einstein_coeff *mol_einst, collisional_transitions *mol_coll, double *mol_popul, int nb_lev, bool acceleration, int verbosity)
{
    bool is_solution_found, is_solution_found_prev;
    int i, lay_nb, nb_cloud_lay, iter_max_nb;
    
    vector<int> bad_layers;
    cloud_layer clayer;

    iteration_control<iteration_scheme_lvg> it_control(it_scheme_lvg);

    if (acceleration)
        iter_max_nb = MAX_NB_ITER_ACC;
    else iter_max_nb = MAX_NB_ITER_EXT;

    nb_cloud_lay = cloud->nb_lay;
    is_solution_found_prev = false;

    for (lay_nb = 0; lay_nb < nb_cloud_lay; lay_nb++) {
        clayer = cloud->lay_array[lay_nb];

        it_scheme_lvg->set_vel_grad(clayer.velg_n);
        it_scheme_lvg->set_dust_parameters(clayer.dust_grain_conc, clayer.dust_grain_temp);
        it_scheme_lvg->set_parameters(clayer.temp_n, clayer.temp_el, clayer.el_conc, clayer.h_conc, clayer.ph2_conc, clayer.oh2_conc, clayer.he_conc,
            clayer.mol_conc, clayer.vel_turb);

        cout << lay_nb << " ";
   
        is_solution_found = false;
        if (lay_nb > 0 && is_solution_found_prev) {
            memcpy(mol_popul + lay_nb * nb_lev, mol_popul + (lay_nb - 1) * nb_lev, nb_lev * sizeof(double));
        }
        else {
            boundary_layer_populations(mol_popul + lay_nb * nb_lev, mol_levels, mol_einst, mol_coll,
                clayer.temp_n, clayer.temp_el, clayer.el_conc, clayer.h_conc, clayer.ph2_conc, clayer.oh2_conc, clayer.he_conc);
        }

        is_solution_found =
            it_control.calculate_populations(mol_popul + lay_nb * nb_lev, iter_max_nb, rel_population_error, acceleration, verbosity);
        
        // trying without acceleration (exception is the methanol)
        if (!is_solution_found && acceleration && mol_levels->mol.name != "CH3OHa" && mol_levels->mol.name != "CH3OHe") {
            acceleration = false;
            iter_max_nb = MAX_NB_ITER_EXT;

            if (lay_nb > 0 && is_solution_found_prev) {
                memcpy(mol_popul + lay_nb * nb_lev, mol_popul + (lay_nb - 1) * nb_lev, nb_lev * sizeof(double));
            }
            else {
                boundary_layer_populations(mol_popul + lay_nb * nb_lev, mol_levels, mol_einst, mol_coll,
                    clayer.temp_n, clayer.temp_el, clayer.el_conc, clayer.h_conc, clayer.ph2_conc, clayer.oh2_conc, clayer.he_conc);
            }

            is_solution_found =
                it_control.calculate_populations(mol_popul + lay_nb * nb_lev, iter_max_nb, rel_population_error, acceleration, verbosity);

            acceleration = true;
            iter_max_nb = MAX_NB_ITER_ACC;
        }

        if (!is_solution_found) {
            bad_layers.push_back(lay_nb);
        }
        is_solution_found_prev = is_solution_found;
    }
    
    cout << endl << "Can not find solution for layers: " << endl;;
    for (i = 0; i < (int)bad_layers.size(); i++) {
        cout << bad_layers[i] << " ";
    }
    cout << endl;
}

void save_mol_data(std::string fname, cloud_data* cloud, double* lev_pop, int nb_lev)
{
	int i, j;
	double a;
	ofstream output;

	output.open(fname.c_str());
    if (!output.is_open())
        cout << "Error in " << SOURCE_NAME << ": can't open file to write molecule population data: " << fname << endl;
    else {
        output << scientific;
        output.precision(4);
        output << left << "! three parameters are given: shock speed (cm/s), turbulent velocity (cm/s), number of specimen levels:" << endl
            << setw(13) << cloud->lay_array[0].vel_n << setw(13) << cloud->lay_array[0].vel_turb << nb_lev << endl;

        output << left << setw(18) << "!depth(cm)" << setw(13) << "gas_temp(K)" << setw(13) << "el_temp" << setw(13) << "dust_temp(K)"
            << setw(13) << "gasvel(cm/s)" << setw(13) << "concHe(cm-3)" << setw(13) << "conc_pH2" << setw(13) << "conc_oH2"
            << setw(13) << "conc_H" << setw(13) << "conc_e" << setw(13) << "conc_H_nucl" << setw(13) << "conc_mol" << "pop_dens";

        for (i = 1; i <= nb_lev; i++) {
            output << left << setw(13) << i;
        }

        for (i = 0; i < cloud->nb_lay; i++) {
            output.precision(10);
            output << left << endl << setw(18) << cloud->lay_array[i].zm;

            output.precision(4);
            output << left << setw(13) << cloud->lay_array[i].temp_n
                << setw(13) << cloud->lay_array[i].temp_el
                << setw(13) << cloud->lay_array[i].av_temp_d
                << setw(13) << cloud->lay_array[i].vel_n
                << setw(13) << cloud->lay_array[i].he_conc
                << setw(13) << cloud->lay_array[i].ph2_conc
                << setw(13) << cloud->lay_array[i].oh2_conc
                << setw(13) << cloud->lay_array[i].h_conc
                << setw(13) << cloud->lay_array[i].el_conc
                << setw(13) << cloud->lay_array[i].tot_h_conc
                << setw(13) << cloud->lay_array[i].mol_conc;

            // population number densities are saving [cm-3]
            for (j = 0; j < nb_lev; j++) {
                a = lev_pop[i * nb_lev + j];
                if (a < 1.e-99)
                    a = 0.;
                output << left << setw(13) << a * cloud->lay_array[i].mol_conc;
            }
        }
        output.close();
    }
}

void save_mol_vibr_data(std::string fname, cloud_data* cloud, energy_diagram* mol_levels, double* lev_pop, int nb_lev)
{
    const int vnb = 6;
    int i, j;
    double a, vdata[vnb];
    ofstream output;

    output.open(fname.c_str());
    if (!output.is_open())
        cout << "Error in " << SOURCE_NAME << ": can't open file to write molecule population data: " << fname << endl;
    else {
        output << scientific;
        output.precision(4);

        output << left << setw(18) << "!depth(cm)" << setw(13) << "gas_temp(K)" << setw(13) << "conc_mol" << "v_dens";
        for (i = 0; i < vnb; i++) {
            output << left << setw(13) << i;
        }

        for (i = 0; i < cloud->nb_lay; i++) {
            output.precision(10);
            output << left << endl << setw(18) << cloud->lay_array[i].zm;

            output.precision(4);
            output << left << setw(13) << cloud->lay_array[i].temp_n << setw(13) << cloud->lay_array[i].mol_conc;

            // vibration state number densities are saved [cm-3]
            memset(vdata, 0, vnb * sizeof(double));
            for (j = 0; j < nb_lev; j++) {
                a = lev_pop[i * nb_lev + j];
                if (mol_levels->lev_array[j].v < vnb)
                    vdata[mol_levels->lev_array[j].v] += a;
            }
            for (j = 0; j < vnb; j++) {
                if (vdata[j] < 1.e-99)
                    vdata[j] = 0.;
                output << left << setw(13) << vdata[j] * cloud->lay_array[i].mol_conc;
            }
        }
        output.close();
    }
}

void calc_oh_masers(string input_data_path, string sim_data_path, string output_path, bool save_line_stat, int verbosity)
{
    bool acceleration, do_simdata_exist;
	int i, isotope, nb_lev_oh, nb_cloud_lay;	
	double mol_mass, spin, dg_ratio, c_abund_pah;
	double * oh_popul;

    stringstream sstr;
    list<transition> trans_list;
    list<transition_data>::const_iterator td_it;

	cloud_layer clayer;
	cloud_data* cloud
		= new cloud_data();
	
    // physical parameters of the particular shock model 
    do_simdata_exist = set_physical_parameters(sim_data_path, cloud);
    do_simdata_exist = set_molecular_conc(sim_data_path, "OH", cloud);
    
    if (do_simdata_exist) {
        // dust model:
        c_abund_pah = 0.;
        dust_model* dust =
            new two_component_dust_model(input_data_path, c_abund_pah, dg_ratio = 0.01, HE_TO_H_NB_RATIO, STANDARD_NB_CR_PHOTONS, verbosity);

        join_layers(cloud, NB_JOINED_LAYERS);
        nb_cloud_lay = cloud->nb_lay;

        cloud->set_vel_turb(MICROTURBULENT_SPEED);
        cloud->set_dust_model(dust);

        // OH with hyperfine levels, 24 or 56
        // up to the HF 20 levels may be involved in the pumping of the maser 1720 MHz (Gray, Proc. IAU Symp. 287, 2012)  
        nb_lev_oh = 56;
        oh_popul = new double[nb_cloud_lay * nb_lev_oh];

        mol_mass = 17. * ATOMIC_MASS_UNIT;
        molecule oh_mol("OH", isotope = 1, mol_mass, spin = 0.5);

        oh_hf_diagram* oh_levels =
            new oh_hf_diagram(input_data_path, oh_mol, nb_lev_oh, verbosity);

        oh_hf_einstein_coeff* oh_einst =
            new oh_hf_einstein_coeff(input_data_path, oh_levels, verbosity);

        oh_hf_collisions* oh_coll =
            new oh_hf_collisions(input_data_path, oh_levels, verbosity);

        iteration_scheme_line_overlap* it_scheme_loverlap
            = new iteration_scheme_line_overlap(input_data_path, verbosity);

        it_scheme_loverlap->init_molecule_data(oh_levels, oh_einst, oh_coll);
        it_scheme_loverlap->set_dust_model(dust);

        calc_molecular_populations(cloud, it_scheme_loverlap, oh_levels, oh_einst, oh_coll, oh_popul, nb_lev_oh, acceleration = false, verbosity);

        transition_data_container* trans_data_oh
            = new transition_data_container(cloud, oh_levels, oh_einst);

        trans_data_oh->set_min_optical_depth(0.01);
        trans_data_oh->find(oh_popul, rel_population_error);

        // adding transitions,
        trans_list.clear();
        oh_trans_list(oh_levels, trans_list);
        trans_data_oh->add(trans_list, oh_popul);

        lim_luminosity_lvg(it_scheme_loverlap, trans_data_oh, cloud, oh_levels, oh_einst, oh_coll, oh_popul);
        trans_data_oh->calc_saturation_depth(BEAMING_FACTOR);

        td_it = trans_data_oh->data.begin();
        while (td_it != trans_data_oh->data.end())
        {
            sstr.clear();
            sstr.str("");
            sstr << output_path + "inv_trans_" + oh_mol.name + "_";
            sstr << (int)(td_it->trans->freq * 1.e-6);
            sstr << ".txt";
            trans_data_oh->save_transition(sstr.str(), *(td_it->trans));
            td_it++;
        }

        sstr.clear();
        sstr.str("");
        sstr << output_path + "inv_trans_" << oh_mol.name << ".txt";
        trans_data_oh->save_full_data(sstr.str());

        sstr.clear();
        sstr.str("");
        sstr << output_path + "pop_dens_" + oh_mol.name + ".txt";
        save_mol_data(sstr.str(), cloud, oh_popul, nb_lev_oh);

        // line statistic for cloud layers in the post shock gas:
        if (save_line_stat) {
            for (i = nb_cloud_lay - 1; i > nb_cloud_lay / 2; i -= 5) {
                clayer = cloud->lay_array[i];
                it_scheme_loverlap->set_vel_grad(clayer.velg_n);
                it_scheme_loverlap->set_dust_parameters(clayer.dust_grain_conc, clayer.dust_grain_temp);
                it_scheme_loverlap->set_parameters(clayer.temp_n, clayer.temp_el, clayer.el_conc, clayer.h_conc, clayer.ph2_conc, clayer.oh2_conc, clayer.he_conc,
                    clayer.mol_conc, clayer.vel_turb);

                sstr.clear();
                sstr.str("");
                sstr << output_path + oh_mol.name + "_line_stat_" << i << ".txt";
                it_scheme_loverlap->calc_line_stat(sstr.str(), oh_popul + i * nb_lev_oh);
            }
        }

        sstr.clear();
        sstr.str("");
        sstr << output_path + "opt_depth_aratio_" + oh_mol.name + ".txt";
        trans_data_oh->save_optical_depth_1(sstr.str());

        sstr.clear();
        sstr.str("");
        sstr << output_path + "opt_depth_freq_" + oh_mol.name + ".txt";
        trans_data_oh->save_optical_depth_2(sstr.str());

        delete[] oh_popul;
        delete oh_levels;
        delete oh_einst;
        delete oh_coll;
        delete trans_data_oh;
        delete it_scheme_loverlap; 
        delete dust;
    }
    delete cloud;
}

void calc_oh_masers_test(string input_data_path, string sim_data_path, string output_path, int verbosity)
{
    bool acceleration, do_simdata_exist;
    int nb_cloud_lay, nb_lev_oh, isotope;
    double c_abund_pah, dg_ratio, mol_mass, spin;
    double* oh_popul;

    stringstream sstr;
    list<transition> trans_list;
    list<transition_data>::const_iterator td_it;

    cloud_data* cloud
        = new cloud_data();    

    // physical parameters of the particular shock model 
    do_simdata_exist = set_physical_parameters(sim_data_path, cloud);
    do_simdata_exist = set_molecular_conc(sim_data_path, "OH", cloud);

    if (do_simdata_exist) {
        // dust model:
        c_abund_pah = 0.;
        dust_model* dust =
            new two_component_dust_model(input_data_path, c_abund_pah, dg_ratio = 0.01, HE_TO_H_NB_RATIO, STANDARD_NB_CR_PHOTONS, verbosity);

        cloud->save_data(output_path + "phys_param1.txt");

        join_layers(cloud, NB_JOINED_LAYERS);
        cloud->save_data(output_path + "phys_param2.txt");

        nb_cloud_lay = cloud->nb_lay;
        cloud->set_vel_turb(MICROTURBULENT_SPEED);
        cloud->set_dust_model(dust);

        // OH with hyperfine levels, 24, 56
        nb_lev_oh = 56;
        oh_popul = new double[nb_cloud_lay * nb_lev_oh];

        mol_mass = 17. * ATOMIC_MASS_UNIT;
        molecule oh_mol("OH", isotope = 1, mol_mass, spin = 0.5);

        oh_hf_diagram* oh_levels =
            new oh_hf_diagram(input_data_path, oh_mol, nb_lev_oh, verbosity);

        oh_hf_einstein_coeff* oh_einst =
            new oh_hf_einstein_coeff(input_data_path, oh_levels, verbosity);

        oh_hf_collisions* oh_coll =
            new oh_hf_collisions(input_data_path, oh_levels, verbosity);

        // no line overlap
        iteration_scheme_lvg* it_scheme_lvg
            = new iteration_scheme_lvg(input_data_path, verbosity);

        it_scheme_lvg->set_dust_model(dust);
        it_scheme_lvg->init_molecule_data(oh_levels, oh_einst, oh_coll);

        calc_molecular_populations(cloud, it_scheme_lvg, oh_levels, oh_einst, oh_coll, oh_popul, nb_lev_oh, acceleration = false, verbosity);

        transition_data_container* trans_data_oh
            = new transition_data_container(cloud, oh_levels, oh_einst);

        trans_data_oh->set_min_optical_depth(0.01);
        trans_data_oh->find(oh_popul, rel_population_error);

        trans_list.clear();
        oh_trans_list(oh_levels, trans_list);
        trans_data_oh->add(trans_list, oh_popul);

        lim_luminosity_lvg(it_scheme_lvg, trans_data_oh, cloud, oh_levels, oh_einst, oh_coll, oh_popul);
        trans_data_oh->calc_saturation_depth(BEAMING_FACTOR);

        td_it = trans_data_oh->data.begin();
        while (td_it != trans_data_oh->data.end())
        {
            sstr.clear();
            sstr.str("");
            sstr << output_path + "inv_trans_" + oh_mol.name + "_";
            sstr << (int)(td_it->trans->freq * 1.e-6);
            sstr << ".txt";
            trans_data_oh->save_transition(sstr.str(), *(td_it->trans));
            td_it++;
        }

        sstr.clear();
        sstr.str("");
        sstr << output_path + "inv_trans_" << oh_mol.name << ".txt";
        trans_data_oh->save_full_data(sstr.str());

        delete[] oh_popul;
        delete oh_levels;
        delete oh_einst;
        delete oh_coll;
        delete trans_data_oh;
        delete it_scheme_lvg;
        delete dust;
    }
    delete cloud;
}

void calc_ch3oh_a_masers(string input_data_path, string sim_data_path, string output_path, bool save_lev_pop, int verbosity)
{
    bool acceleration, do_simdata_exist;  
    int isotope, nb_lev_ch3oh, nb_vibr_ch3oh, ang_mom_max, nb_cloud_lay;
    double dg_ratio, c_abund_pah, mol_mass, spin;
    double* ch3oh_a_popul;

    stringstream sstr;
    list<transition> trans_list;
    list<transition_data>::const_iterator td_it;

    cloud_data* cloud
        = new cloud_data();

    // physical parameters of the particular shock model 
    do_simdata_exist = set_physical_parameters(sim_data_path, cloud);
    do_simdata_exist = set_molecular_conc(sim_data_path, "CH3OH", cloud, a_meth_fraction);
         
    if (do_simdata_exist) {
        // dust model:
        c_abund_pah = 0.;
        const dust_model* dust = 
            new two_component_dust_model(input_data_path, c_abund_pah, dg_ratio = 0.01, HE_TO_H_NB_RATIO, STANDARD_NB_CR_PHOTONS, verbosity);

        join_layers(cloud, NB_JOINED_LAYERS);
        nb_cloud_lay = cloud->nb_lay;

        cloud->set_vel_turb(MICROTURBULENT_SPEED);
        cloud->set_dust_model(dust);
    
        // The calculation of the A-class methanol populations;
        // Max nb of levels for states vt=0,1,2 for which spectroscopic data are available is 1455 (angular momentum <= 22); 
        // collisional coefficients are available for 3*256 levels belonging to vt=0,1,2 and having angular momentum <=15;
        nb_lev_ch3oh = 768;
        nb_vibr_ch3oh = 2;
        ang_mom_max = 15;
        ch3oh_a_popul = new double[nb_cloud_lay * nb_lev_ch3oh];

        mol_mass = 32. * ATOMIC_MASS_UNIT;
        molecule ch3oh_a_mol("CH3OHa", isotope = 1, mol_mass, spin = 1.5);

        ch3oh_diagram* ch3oh_a_levels =
            new ch3oh_diagram(input_data_path, ch3oh_a_mol, nb_lev_ch3oh, nb_vibr_ch3oh, ang_mom_max, verbosity);

        ch3oh_einstein_coeff* ch3oh_a_einst
            = new ch3oh_einstein_coeff(input_data_path, ch3oh_a_levels, verbosity);

        ch3oh_collisions* ch3oh_a_coll
            = new ch3oh_collisions(input_data_path, ch3oh_a_levels, verbosity);

        iteration_scheme_lvg* it_scheme_lvg
            = new iteration_scheme_lvg(input_data_path, verbosity);

        it_scheme_lvg->set_dust_model(dust);
        it_scheme_lvg->init_molecule_data(ch3oh_a_levels, ch3oh_a_einst, ch3oh_a_coll);

        calc_molecular_populations(cloud, it_scheme_lvg, ch3oh_a_levels, ch3oh_a_einst, ch3oh_a_coll, ch3oh_a_popul, nb_lev_ch3oh, acceleration = true, verbosity);

        transition_data_container* trans_data_ch3oh_a
            = new transition_data_container(cloud, ch3oh_a_levels, ch3oh_a_einst);
        
        trans_data_ch3oh_a->set_min_optical_depth(0.1);
        trans_data_ch3oh_a->find(ch3oh_a_popul, rel_population_error);
        
        // the transitions in the given list; 
        trans_list.clear();
        ch3oh_classI_trans_list(ch3oh_a_levels, trans_list);  // only class I masers
        trans_data_ch3oh_a->add(trans_list, ch3oh_a_popul);

        lim_luminosity_lvg(it_scheme_lvg, trans_data_ch3oh_a, cloud, ch3oh_a_levels, ch3oh_a_einst, ch3oh_a_coll, ch3oh_a_popul);
        trans_data_ch3oh_a->calc_saturation_depth(BEAMING_FACTOR);

        // the transitions that have population inversion and are in the given list, are saved individually,
        td_it = trans_data_ch3oh_a->data.begin();
        while (td_it != trans_data_ch3oh_a->data.end())
        {
            sstr.clear();
            sstr.str("");
            sstr << output_path + "inv_trans_" + ch3oh_a_mol.name + "_";
            sstr << (int)(td_it->trans->freq * 1.e-6);
            sstr << ".txt";
            trans_data_ch3oh_a->save_transition(sstr.str(), *(td_it->trans));
            td_it++;
        }
        
        sstr.clear();
        sstr.str("");
        sstr << output_path + "inv_trans_" << ch3oh_a_mol.name << ".txt";
        trans_data_ch3oh_a->save_full_data(sstr.str());

        if (save_lev_pop) {
            sstr.clear();
            sstr.str("");
            sstr << output_path + "pop_dens_" + ch3oh_a_mol.name + ".txt";
            save_mol_data(sstr.str(), cloud, ch3oh_a_popul, nb_lev_ch3oh);

            sstr.clear();
            sstr.str("");
            sstr << output_path + "vibr_dens_" + ch3oh_a_mol.name + ".txt";
            save_mol_vibr_data(sstr.str(), cloud, ch3oh_a_levels, ch3oh_a_popul, nb_lev_ch3oh);
        }

        sstr.clear();
        sstr.str("");
        sstr << output_path + "opt_depth_aratio_" + ch3oh_a_mol.name + ".txt";
        trans_data_ch3oh_a->save_optical_depth_1(sstr.str());

        sstr.clear();
        sstr.str("");
        sstr << output_path + "opt_depth_freq_" + ch3oh_a_mol.name + ".txt";
        trans_data_ch3oh_a->save_optical_depth_2(sstr.str());

        delete[] ch3oh_a_popul;
        delete ch3oh_a_levels;
        delete ch3oh_a_einst;
        delete ch3oh_a_coll;
        delete trans_data_ch3oh_a;
        delete it_scheme_lvg; 
        delete dust;
    }
    delete cloud;
}

void calc_ch3oh_e_masers(string input_data_path, string sim_data_path, string output_path, int verbosity)
{
    bool acceleration, do_simdata_exist;
    int isotope, nb_lev_ch3oh, nb_vibr_ch3oh, ang_mom_max, nb_cloud_lay;
    double dg_ratio, c_abund_pah, mol_mass, spin;
    double * ch3oh_e_popul;

    stringstream sstr;
    list<transition> trans_list;
    list<transition_data>::const_iterator td_it;

    cloud_data* cloud
        = new cloud_data();

    // physical parameters of the particular shock model 
    do_simdata_exist = set_physical_parameters(sim_data_path, cloud);
    do_simdata_exist = set_molecular_conc(sim_data_path, "CH3OH", cloud, 1. - a_meth_fraction);
        
    if (do_simdata_exist) {
        // dust model:
        c_abund_pah = 0.;
        const dust_model* dust =
            new two_component_dust_model(input_data_path, c_abund_pah, dg_ratio = 0.01, HE_TO_H_NB_RATIO, STANDARD_NB_CR_PHOTONS, verbosity);

        join_layers(cloud, NB_JOINED_LAYERS);
        nb_cloud_lay = cloud->nb_lay;

        cloud->set_vel_turb(MICROTURBULENT_SPEED);
        cloud->set_dust_model(dust);
               
        // The calculation of the E-class methanol populations;
        nb_lev_ch3oh = 768;
        nb_vibr_ch3oh = 2;
        ang_mom_max = 15;
        ch3oh_e_popul = new double[nb_cloud_lay * nb_lev_ch3oh];

        mol_mass = 32. * ATOMIC_MASS_UNIT;
        molecule ch3oh_e_mol("CH3OHe", isotope = 1, mol_mass, spin = 0.5);

        ch3oh_diagram* ch3oh_e_levels =
            new ch3oh_diagram(input_data_path, ch3oh_e_mol, nb_lev_ch3oh, nb_vibr_ch3oh, ang_mom_max, verbosity);

        ch3oh_einstein_coeff* ch3oh_e_einst
            = new ch3oh_einstein_coeff(input_data_path, ch3oh_e_levels, verbosity);

        ch3oh_collisions* ch3oh_e_coll
            = new ch3oh_collisions(input_data_path, ch3oh_e_levels, verbosity);

                    
        iteration_scheme_lvg* it_scheme_lvg
            = new iteration_scheme_lvg(input_data_path, verbosity);

        it_scheme_lvg->set_dust_model(dust);
        it_scheme_lvg->init_molecule_data(ch3oh_e_levels, ch3oh_e_einst, ch3oh_e_coll);

        calc_molecular_populations(cloud, it_scheme_lvg, ch3oh_e_levels, ch3oh_e_einst, ch3oh_e_coll, ch3oh_e_popul, nb_lev_ch3oh, acceleration = true, verbosity);

        transition_data_container* trans_data_ch3oh_e
            = new transition_data_container(cloud, ch3oh_e_levels, ch3oh_e_einst);

        trans_data_ch3oh_e->set_min_optical_depth(0.1);
        trans_data_ch3oh_e->find(ch3oh_e_popul, rel_population_error);
        
        trans_list.clear();
        ch3oh_classI_trans_list(ch3oh_e_levels, trans_list);
        trans_data_ch3oh_e->add(trans_list, ch3oh_e_popul);

        lim_luminosity_lvg(it_scheme_lvg, trans_data_ch3oh_e, cloud, ch3oh_e_levels, ch3oh_e_einst, ch3oh_e_coll, ch3oh_e_popul);
        trans_data_ch3oh_e->calc_saturation_depth(BEAMING_FACTOR);

        td_it = trans_data_ch3oh_e->data.begin();
        while (td_it != trans_data_ch3oh_e->data.end())
        {
            sstr.clear();
            sstr.str("");
            sstr << output_path + "inv_trans_" + ch3oh_e_mol.name + "_";
            sstr << (int)(td_it->trans->freq * 1.e-6);
            sstr << ".txt";
            trans_data_ch3oh_e->save_transition(sstr.str(), *(td_it->trans));  
            td_it++;
        }
        
        sstr.clear();
        sstr.str("");
        sstr << output_path + "inv_trans_" << ch3oh_e_mol.name << ".txt";
        trans_data_ch3oh_e->save_full_data(sstr.str());

        sstr.clear();
        sstr.str("");
        sstr << output_path + "opt_depth_aratio_" + ch3oh_e_mol.name + ".txt";
        trans_data_ch3oh_e->save_optical_depth_1(sstr.str());

        sstr.clear();
        sstr.str("");
        sstr << output_path + "opt_depth_freq_" + ch3oh_e_mol.name + ".txt";
        trans_data_ch3oh_e->save_optical_depth_2(sstr.str());

        delete[] ch3oh_e_popul;
        delete ch3oh_e_levels;
        delete ch3oh_e_einst;
        delete ch3oh_e_coll;
        delete trans_data_ch3oh_e;
        delete it_scheme_lvg; 
        delete dust;
    }
    delete cloud;
}

// check acceleration mode and nb of iterations, 
void calc_h2o_para_masers(string input_data_path, string sim_data_path, string output_path, bool save_lev_pop, bool is_shock_data_used, int verbosity)
{
    bool acceleration, do_simdata_exist;
    int isotope, nb_lev_h2o, nb_vibr_h2o, nb_cloud_lay;
    double dg_ratio, c_abund_pah, mol_mass, spin;
    double* ph2o_popul;

    stringstream sstr;
    list<transition> trans_list;
    list<transition_data>::const_iterator td_it;

    cloud_data* cloud
        = new cloud_data();

    // physical parameters of the particular shock model 
    do_simdata_exist = set_physical_parameters(sim_data_path, cloud);
    do_simdata_exist = set_molecular_conc(sim_data_path, "H2O", cloud, 1. - oh2o_fraction);

    if (do_simdata_exist) {
        // dust model:
        c_abund_pah = 0.;
        const dust_model* dust =
            new two_component_dust_model(input_data_path, c_abund_pah, dg_ratio = 0.01, HE_TO_H_NB_RATIO, STANDARD_NB_CR_PHOTONS, verbosity);

        join_layers(cloud, NB_JOINED_LAYERS);
        nb_cloud_lay = cloud->nb_lay;

        cloud->set_vel_turb(MICROTURBULENT_SPEED);
        cloud->set_dust_model(dust);
    
        // The calculation of the para-H2O populations;
        // H2O molecule data:
        nb_lev_h2o = 150;
        nb_vibr_h2o = 1;
        ph2o_popul = new double[nb_cloud_lay * nb_lev_h2o];

        mol_mass = 18. * ATOMIC_MASS_UNIT;
        molecule ph2o_mol("pH2O", isotope = 1, mol_mass, spin = 0.);
        
        h2o_diagram* ph2o_levels 
            = new h2o_diagram(input_data_path, ph2o_mol, nb_lev_h2o, nb_vibr_h2o, verbosity);
       
        // the pointer to derivative class must be used here (not to base class):
        h2o_einstein_coeff *ph2o_einst 
            = new h2o_einstein_coeff(input_data_path, ph2o_levels, verbosity);
        
        h2o_collisions *ph2o_coll 
            = new h2o_collisions(input_data_path, ph2o_levels, false, verbosity);

        iteration_scheme_lvg* it_scheme_lvg
                = new iteration_scheme_lvg(input_data_path, verbosity);

        it_scheme_lvg->set_dust_model(dust);
        it_scheme_lvg->init_molecule_data(ph2o_levels, ph2o_einst, ph2o_coll);

        if (!is_shock_data_used) {
            calc_molecular_populations(cloud, it_scheme_lvg, ph2o_levels, ph2o_einst, ph2o_coll, ph2o_popul, nb_lev_h2o, acceleration = true, verbosity);
        }
        else {
            set_level_pop(sim_data_path + "sim_data_ph2o.txt", cloud, nb_lev_h2o, ph2o_popul);
        }

        transition_data_container* trans_data_ph2o
            = new transition_data_container(cloud, ph2o_levels, ph2o_einst);

        trans_data_ph2o->set_min_optical_depth(0.1);
        trans_data_ph2o->find(ph2o_popul, rel_population_error);
        
        trans_list.clear();
        para_h2o_trans_list(ph2o_levels, trans_list);
        trans_data_ph2o->add(trans_list, ph2o_popul);

        lim_luminosity_lvg(it_scheme_lvg, trans_data_ph2o, cloud, ph2o_levels, ph2o_einst, ph2o_coll, ph2o_popul);
        trans_data_ph2o->calc_saturation_depth(BEAMING_FACTOR);

        td_it = trans_data_ph2o->data.begin();
        while (td_it != trans_data_ph2o->data.end())
        {
            sstr.clear();
            sstr.str("");
            sstr << output_path + "inv_trans_" + ph2o_mol.name + "_";
            sstr << (int)(td_it->trans->freq * 1.e-6);
            sstr << ".txt";
            trans_data_ph2o->save_transition(sstr.str(), *(td_it->trans));
            td_it++;
        }
        
        sstr.clear();
        sstr.str("");
        sstr << output_path + "inv_trans_" << ph2o_mol.name << ".txt";
        trans_data_ph2o->save_full_data(sstr.str());

        if (save_lev_pop) {
            sstr.clear();
            sstr.str("");
            sstr << output_path + "pop_dens_" + ph2o_mol.name + ".txt";
            save_mol_data(sstr.str(), cloud, ph2o_popul, nb_lev_h2o);
        }

        delete[] ph2o_popul;
        delete ph2o_levels;
        delete ph2o_einst;
        delete ph2o_coll;
        delete trans_data_ph2o;
        delete it_scheme_lvg;
        delete dust;
    }
    delete cloud;
}

// check acceleration mode and nb of iterations, 
void calc_h2o_ortho_masers(string input_data_path, string sim_data_path, string output_path, bool save_lev_pop, int verbosity)
{
    bool acceleration, do_simdata_exist;
    int isotope, nb_lev_h2o, nb_vibr_h2o, nb_cloud_lay;
    double dg_ratio, c_abund_pah, mol_mass, spin;
    double* oh2o_popul;

    stringstream sstr;
    list<transition> trans_list;
    list<transition_data>::const_iterator td_it;

    cloud_data* cloud
        = new cloud_data();

    // physical parameters of the particular shock model 
    do_simdata_exist = set_physical_parameters(sim_data_path, cloud);
    do_simdata_exist = set_molecular_conc(sim_data_path, "H2O", cloud, oh2o_fraction);

    if (do_simdata_exist) {
        // dust model:
        c_abund_pah = 0.;
        const dust_model* dust =
            new two_component_dust_model(input_data_path, c_abund_pah, dg_ratio = 0.01, HE_TO_H_NB_RATIO, STANDARD_NB_CR_PHOTONS, verbosity);

        join_layers(cloud, NB_JOINED_LAYERS);
        nb_cloud_lay = cloud->nb_lay;

        cloud->set_vel_turb(MICROTURBULENT_SPEED);
        cloud->set_dust_model(dust);

        // The calculation of the ortho-H2O populations;
        // H2O molecule data:
        nb_lev_h2o = 150;
        nb_vibr_h2o = 1;
        oh2o_popul = new double[nb_cloud_lay * nb_lev_h2o];

        mol_mass = 18. * ATOMIC_MASS_UNIT;
        molecule oh2o_mol("oH2O", isotope = 1, mol_mass, spin = 1.);

        h2o_diagram* oh2o_levels
            = new h2o_diagram(input_data_path, oh2o_mol, nb_lev_h2o, nb_vibr_h2o, verbosity);

        // the pointer to derivative class must be used here (not to base class):
        h2o_einstein_coeff* oh2o_einst
            = new h2o_einstein_coeff(input_data_path, oh2o_levels, verbosity);

        h2o_collisions* oh2o_coll
            = new h2o_collisions(input_data_path, oh2o_levels, false, verbosity);

        iteration_scheme_lvg* it_scheme_lvg
            = new iteration_scheme_lvg(input_data_path, verbosity);

        it_scheme_lvg->set_dust_model(dust);
        it_scheme_lvg->init_molecule_data(oh2o_levels, oh2o_einst, oh2o_coll);

        calc_molecular_populations(cloud, it_scheme_lvg, oh2o_levels, oh2o_einst, oh2o_coll, oh2o_popul, nb_lev_h2o, acceleration = true, verbosity);
       
        transition_data_container* trans_data_oh2o
            = new transition_data_container(cloud, oh2o_levels, oh2o_einst);

        trans_data_oh2o->set_min_optical_depth(0.1);
        trans_data_oh2o->find(oh2o_popul, rel_population_error);
        
        trans_list.clear();
        ortho_h2o_trans_list(oh2o_levels, trans_list);
        trans_data_oh2o->add(trans_list, oh2o_popul);
        
        lim_luminosity_lvg(it_scheme_lvg, trans_data_oh2o, cloud, oh2o_levels, oh2o_einst, oh2o_coll, oh2o_popul);
        trans_data_oh2o->calc_saturation_depth(BEAMING_FACTOR);

        td_it = trans_data_oh2o->data.begin();
        while (td_it != trans_data_oh2o->data.end())
        {
            sstr.clear();
            sstr.str("");
            sstr << output_path + "inv_trans_" + oh2o_mol.name + "_";
            sstr << (int)(td_it->trans->freq * 1.e-6);
            sstr << ".txt";
            trans_data_oh2o->save_transition(sstr.str(), *(td_it->trans));
            td_it++;
        }       

        sstr.clear();
        sstr.str("");
        sstr << output_path + "inv_trans_" << oh2o_mol.name << ".txt";
        trans_data_oh2o->save_full_data(sstr.str());

        if (save_lev_pop) {
            sstr.clear();
            sstr.str("");
            sstr << output_path + "pop_dens_" + oh2o_mol.name + ".txt";
            save_mol_data(sstr.str(), cloud, oh2o_popul, nb_lev_h2o);
        }

        delete[] oh2o_popul;
        delete oh2o_levels;
        delete oh2o_einst;
        delete oh2o_coll;
        delete trans_data_oh2o;
        delete it_scheme_lvg;
        delete dust;
    }
    delete cloud;
}

void calc_nh3_ortho_masers(string input_data_path, string sim_data_path, string output_path, int verbosity)
{
    bool acceleration, do_simdata_exist;
    int isotope, nb_lev_onh3, nb_cloud_lay;
    double dg_ratio, c_abund_pah, mol_mass, spin;
    double* onh3_popul;

    stringstream sstr;
    list<transition> trans_list;
    list<transition_data>::const_iterator td_it;

    cloud_data* cloud
        = new cloud_data();

    // physical parameters of the particular shock model 
    do_simdata_exist = set_physical_parameters(sim_data_path, cloud);
    do_simdata_exist = set_molecular_conc(sim_data_path, "NH3", cloud, onh3_fraction);

    if (do_simdata_exist) {
        // dust model:
        c_abund_pah = 0.;
        const dust_model* dust =
            new two_component_dust_model(input_data_path, c_abund_pah, dg_ratio = 0.01, HE_TO_H_NB_RATIO, STANDARD_NB_CR_PHOTONS, verbosity);

        join_layers(cloud, NB_JOINED_LAYERS);
        nb_cloud_lay = cloud->nb_lay;

        cloud->set_vel_turb(MICROTURBULENT_SPEED);
        cloud->set_dust_model(dust);

        // ortho-NH3
        // NH3 molecule data, o-NH3 has k = 3n, n is an integer, for p-NH3 k != 3n
        nb_lev_onh3 = 22; // ortho-NH3: He coll data - 22, H2 coll data - 17
        onh3_popul = new double[nb_cloud_lay * nb_lev_onh3];

        mol_mass = 17. * ATOMIC_MASS_UNIT;
        molecule onh3_mol("oNH3", isotope = 1, mol_mass, spin = 1.5);

        nh3_diagram* onh3_levels =
            new nh3_diagram(input_data_path, onh3_mol, nb_lev_onh3, verbosity);

        nh3_einstein_coeff* onh3_einst =
            new nh3_einstein_coeff(input_data_path, onh3_levels, verbosity);

        nh3_collisions* onh3_coll =
            new nh3_collisions(input_data_path, onh3_levels, verbosity);

        iteration_scheme_lvg* it_scheme_lvg
            = new iteration_scheme_lvg(input_data_path, verbosity);

        it_scheme_lvg->set_dust_model(dust);
        it_scheme_lvg->init_molecule_data(onh3_levels, onh3_einst, onh3_coll);

        calc_molecular_populations(cloud, it_scheme_lvg, onh3_levels, onh3_einst, onh3_coll, onh3_popul, nb_lev_onh3, acceleration = true, verbosity);

        transition_data_container* trans_data_onh3
            = new transition_data_container(cloud, onh3_levels, onh3_einst);

        trans_data_onh3->set_min_optical_depth(0.1);
        trans_data_onh3->find(onh3_popul, rel_population_error);
        
        trans_list.clear();
        nh3_trans_list(onh3_levels, trans_list);
        trans_data_onh3->add(trans_list, onh3_popul);

        lim_luminosity_lvg(it_scheme_lvg, trans_data_onh3, cloud, onh3_levels, onh3_einst, onh3_coll, onh3_popul);
        trans_data_onh3->calc_saturation_depth(BEAMING_FACTOR);

        td_it = trans_data_onh3->data.begin();
        while (td_it != trans_data_onh3->data.end())
        {
            sstr.clear();
            sstr.str("");
            sstr << output_path + "inv_trans_" + onh3_mol.name + "_";
            sstr << (int)(td_it->trans->freq * 1.e-6);
            sstr << ".txt";
            trans_data_onh3->save_transition(sstr.str(), *(td_it->trans));
            td_it++;
        }

        sstr.clear();
        sstr.str("");
        sstr << output_path + "inv_trans_" << onh3_mol.name << ".txt";
        trans_data_onh3->save_full_data(sstr.str());

        delete[] onh3_popul;
        delete onh3_levels;
        delete onh3_einst;
        delete onh3_coll;
        delete trans_data_onh3;
        delete it_scheme_lvg;
        delete dust;
    }
    delete cloud;
}

void calc_nh3_para_masers(string input_data_path, string sim_data_path, string output_path, int verbosity)
{
    bool acceleration, do_simdata_exist;
    int isotope, nb_lev_pnh3, nb_cloud_lay;
    double dg_ratio, c_abund_pah, mol_mass, spin;
    double* pnh3_popul;

    stringstream sstr;
    list<transition> trans_list;
    list<transition_data>::const_iterator td_it;

    cloud_data* cloud
        = new cloud_data();

    // physical parameters of the particular shock model 
    do_simdata_exist = set_physical_parameters(sim_data_path, cloud);
    do_simdata_exist = set_molecular_conc(sim_data_path, "NH3", cloud, 1. - onh3_fraction);

    if (do_simdata_exist) {
        // dust model:
        c_abund_pah = 0.;
        const dust_model* dust =
            new two_component_dust_model(input_data_path, c_abund_pah, dg_ratio = 0.01, HE_TO_H_NB_RATIO, STANDARD_NB_CR_PHOTONS, verbosity);

        join_layers(cloud, NB_JOINED_LAYERS);
        nb_cloud_lay = cloud->nb_lay;

        cloud->set_vel_turb(MICROTURBULENT_SPEED);
        cloud->set_dust_model(dust);

        // para-NH3
        nb_lev_pnh3 = 34; // para-NH3: He coll data - 16, H2 coll data - 34
        pnh3_popul = new double[nb_cloud_lay * nb_lev_pnh3];

        mol_mass = 17. * ATOMIC_MASS_UNIT;
        molecule pnh3_mol("pNH3", isotope = 1, mol_mass, spin = 0.5);

        nh3_diagram* pnh3_levels =
            new nh3_diagram(input_data_path, pnh3_mol, nb_lev_pnh3, verbosity);

        nh3_einstein_coeff* pnh3_einst =
            new nh3_einstein_coeff(input_data_path, pnh3_levels, verbosity);

        nh3_collisions* pnh3_coll =
            new nh3_collisions(input_data_path, pnh3_levels, verbosity);


        iteration_scheme_lvg* it_scheme_lvg
            = new iteration_scheme_lvg(input_data_path, verbosity);

        it_scheme_lvg->set_dust_model(dust);
        it_scheme_lvg->init_molecule_data(pnh3_levels, pnh3_einst, pnh3_coll);

        calc_molecular_populations(cloud, it_scheme_lvg, pnh3_levels, pnh3_einst, pnh3_coll, pnh3_popul, nb_lev_pnh3, acceleration = true, verbosity);

        transition_data_container* trans_data_pnh3
            = new transition_data_container(cloud, pnh3_levels, pnh3_einst);

        trans_data_pnh3->set_min_optical_depth(0.1);
        trans_data_pnh3->find(pnh3_popul, rel_population_error);
  
        trans_list.clear();
        nh3_trans_list(pnh3_levels, trans_list);
        trans_data_pnh3->add(trans_list, pnh3_popul);

        lim_luminosity_lvg(it_scheme_lvg, trans_data_pnh3, cloud, pnh3_levels, pnh3_einst, pnh3_coll, pnh3_popul);
        trans_data_pnh3->calc_saturation_depth(BEAMING_FACTOR);

        td_it = trans_data_pnh3->data.begin();
        while (td_it != trans_data_pnh3->data.end())
        {
            sstr.clear();
            sstr.str("");
            sstr << output_path + "inv_trans_" + pnh3_mol.name + "_";
            sstr << (int)(td_it->trans->freq * 1.e-6);
            sstr << ".txt";
            trans_data_pnh3->save_transition(sstr.str(), *(td_it->trans));
            td_it++;
        }

        sstr.clear();
        sstr.str("");
        sstr << output_path + "inv_trans_" << pnh3_mol.name << ".txt";
        trans_data_pnh3->save_full_data(sstr.str());

        delete[] pnh3_popul;
        delete pnh3_levels;
        delete pnh3_einst;
        delete pnh3_coll;
        delete trans_data_pnh3;
        delete it_scheme_lvg;
        delete dust;
    }
    delete cloud;
}

void calc_h2co_ortho_masers(string input_data_path, string sim_data_path, string output_path, int verbosity)
{
    bool acceleration, do_simdata_exist;
    int isotope, nb_lev_oh2co, nb_cloud_lay;
    double dg_ratio, c_abund_pah, mol_mass, spin;
    double* oh2co_popul;

    stringstream sstr;
    list<transition> trans_list;
    list<transition_data>::const_iterator td_it;

    cloud_data* cloud
        = new cloud_data();

    // physical parameters of the particular shock model 
    do_simdata_exist = set_physical_parameters(sim_data_path, cloud);
    do_simdata_exist = set_molecular_conc(sim_data_path, "H2CO", cloud, oh2co_fraction);

    if (do_simdata_exist) {
        // dust model:
        c_abund_pah = 0.;
        const dust_model* dust =
            new two_component_dust_model(input_data_path, c_abund_pah, dg_ratio = 0.01, HE_TO_H_NB_RATIO, STANDARD_NB_CR_PHOTONS, verbosity);

        join_layers(cloud, NB_JOINED_LAYERS);
        nb_cloud_lay = cloud->nb_lay;

        cloud->set_vel_turb(MICROTURBULENT_SPEED);
        cloud->set_dust_model(dust);

        // ortho-H2CO
        nb_lev_oh2co = 40; // ortho-H2CO: H2 coll data - 40
        oh2co_popul = new double[nb_cloud_lay * nb_lev_oh2co];

        mol_mass = 30. * ATOMIC_MASS_UNIT;
        molecule oh2co_mol("oH2CO", isotope = 1, mol_mass, spin = 1);  // ortho-H2CO, spin = 1

        h2co_diagram* oh2co_levels =
            new h2co_diagram(input_data_path, oh2co_mol, nb_lev_oh2co, verbosity);

        h2co_einstein_coeff* oh2co_einst =
            new h2co_einstein_coeff(input_data_path, oh2co_levels, verbosity);
        
        h2co_collisions* oh2co_coll =
            new h2co_collisions(input_data_path, oh2co_levels, verbosity);

        iteration_scheme_lvg* it_scheme_lvg
            = new iteration_scheme_lvg(input_data_path, verbosity);

        it_scheme_lvg->set_dust_model(dust);
        it_scheme_lvg->init_molecule_data(oh2co_levels, oh2co_einst, oh2co_coll);

        calc_molecular_populations(cloud, it_scheme_lvg, oh2co_levels, oh2co_einst, oh2co_coll, oh2co_popul, nb_lev_oh2co, acceleration = true, verbosity);

        transition_data_container* trans_data_oh2co
            = new transition_data_container(cloud, oh2co_levels, oh2co_einst);

        trans_data_oh2co->set_min_optical_depth(0.1);
        trans_data_oh2co->find(oh2co_popul, rel_population_error);

        trans_list.clear();
        h2co_trans_list(oh2co_levels, trans_list);
        trans_data_oh2co->add(trans_list, oh2co_popul);

        lim_luminosity_lvg(it_scheme_lvg, trans_data_oh2co, cloud, oh2co_levels, oh2co_einst, oh2co_coll, oh2co_popul);
        trans_data_oh2co->calc_saturation_depth(BEAMING_FACTOR);

        td_it = trans_data_oh2co->data.begin();
        while (td_it != trans_data_oh2co->data.end())
        {
            sstr.clear();
            sstr.str("");
            sstr << output_path + "inv_trans_" + oh2co_mol.name + "_";
            sstr << (int)(td_it->trans->freq * 1.e-6);
            sstr << ".txt";
            trans_data_oh2co->save_transition(sstr.str(), *(td_it->trans));
            td_it++;
        }

        sstr.clear();
        sstr.str("");
        sstr << output_path + "inv_trans_" << oh2co_mol.name << ".txt";
        trans_data_oh2co->save_full_data(sstr.str());

        delete[] oh2co_popul;
        delete oh2co_levels;
        delete oh2co_einst;
        delete oh2co_coll;
        delete trans_data_oh2co;
        delete it_scheme_lvg;
        delete dust;
    }
    delete cloud;
}

void calc_h2co_para_masers(string input_data_path, string sim_data_path, string output_path, int verbosity)
{
    bool acceleration, do_simdata_exist;
    int isotope, nb_lev_ph2co, nb_cloud_lay;
    double dg_ratio, c_abund_pah, mol_mass, spin;
    double* ph2co_popul;

    stringstream sstr;
    list<transition> trans_list;
    list<transition_data>::const_iterator td_it;

    cloud_data* cloud
        = new cloud_data();

    // physical parameters of the particular shock model 
    do_simdata_exist = set_physical_parameters(sim_data_path, cloud);
    do_simdata_exist = set_molecular_conc(sim_data_path, "H2CO", cloud, 1. - oh2co_fraction);

    if (do_simdata_exist) {
        // dust model:
        c_abund_pah = 0.;
        const dust_model* dust =
            new two_component_dust_model(input_data_path, c_abund_pah, dg_ratio = 0.01, HE_TO_H_NB_RATIO, STANDARD_NB_CR_PHOTONS, verbosity);

        join_layers(cloud, NB_JOINED_LAYERS);
        nb_cloud_lay = cloud->nb_lay;

        cloud->set_vel_turb(MICROTURBULENT_SPEED);
        cloud->set_dust_model(dust);

        // para-H2CO
        nb_lev_ph2co = 41; // para-H2CO: H2 coll data - 41
        ph2co_popul = new double[nb_cloud_lay * nb_lev_ph2co];

        mol_mass = 30. * ATOMIC_MASS_UNIT;
        molecule ph2co_mol("pH2CO", isotope = 1, mol_mass, spin = 0);  // para-H2CO, spin = 0

        h2co_diagram* ph2co_levels =
            new h2co_diagram(input_data_path, ph2co_mol, nb_lev_ph2co, verbosity);

        h2co_einstein_coeff* ph2co_einst =
            new h2co_einstein_coeff(input_data_path, ph2co_levels, verbosity);

        h2co_collisions* ph2co_coll =
            new h2co_collisions(input_data_path, ph2co_levels, verbosity);

        iteration_scheme_lvg* it_scheme_lvg
            = new iteration_scheme_lvg(input_data_path, verbosity);

        it_scheme_lvg->set_dust_model(dust);
        it_scheme_lvg->init_molecule_data(ph2co_levels, ph2co_einst, ph2co_coll);

        calc_molecular_populations(cloud, it_scheme_lvg, ph2co_levels, ph2co_einst, ph2co_coll, ph2co_popul, nb_lev_ph2co, acceleration = true, verbosity);

        transition_data_container* trans_data_ph2co
            = new transition_data_container(cloud, ph2co_levels, ph2co_einst);

        trans_data_ph2co->set_min_optical_depth(0.1);
        trans_data_ph2co->find(ph2co_popul, rel_population_error);

        trans_list.clear();
        h2co_trans_list(ph2co_levels, trans_list);
        trans_data_ph2co->add(trans_list, ph2co_popul);

        lim_luminosity_lvg(it_scheme_lvg, trans_data_ph2co, cloud, ph2co_levels, ph2co_einst, ph2co_coll, ph2co_popul);
        trans_data_ph2co->calc_saturation_depth(BEAMING_FACTOR);

        td_it = trans_data_ph2co->data.begin();
        while (td_it != trans_data_ph2co->data.end())
        {
            sstr.clear();
            sstr.str("");
            sstr << output_path + "inv_trans_" + ph2co_mol.name + "_";
            sstr << (int)(td_it->trans->freq * 1.e-6);
            sstr << ".txt";
            trans_data_ph2co->save_transition(sstr.str(), *(td_it->trans));
            td_it++;
        }

        sstr.clear();
        sstr.str("");
        sstr << output_path + "inv_trans_" << ph2co_mol.name << ".txt";
        trans_data_ph2co->save_full_data(sstr.str());

        delete[] ph2co_popul;
        delete ph2co_levels;
        delete ph2co_einst;
        delete ph2co_coll;
        delete trans_data_ph2co;
        delete it_scheme_lvg;
        delete dust;
    }
    delete cloud;
}

void calc_nb_line_overlaps(string input_data_path, string output_path)
{
    int verbosity = 1;
    int isotope, nb_lev_oh, nb_d, nb_tr;
    double mol_mass, spin, vel_width;

    string fname;
    ofstream outfile;

    nb_lev_oh = 24;
    mol_mass = 17. * ATOMIC_MASS_UNIT;

    molecule oh_mol("OH", isotope = 1, mol_mass, spin = 0.5);

    oh_hf_diagram* oh_levels =
        new oh_hf_diagram(input_data_path, oh_mol, nb_lev_oh, verbosity);

    oh_hf_einstein_coeff* oh_einst =
        new oh_hf_einstein_coeff(input_data_path, oh_levels, verbosity);

    fname = output_path + "OH_overlap.txt";
    outfile.open(fname.c_str());
    outfile.setf(ios::scientific);
    outfile << setprecision(4);

    outfile << "velocity width (cm/s) - double overlap - triple overlap" << endl;
    for (vel_width = 1.e+4; vel_width < 1.e+5; vel_width *= 1.05) {
        get_nb_overlap_lines(oh_levels, oh_einst, vel_width, nb_d, nb_tr);
        outfile << left << setw(14) << vel_width << setw(10) << nb_d << setw(10) << nb_tr << endl;
    }
    outfile.close();
}

void ch3oh_coll_extrapolations(string input_data_path, string output_path)
{
    int verbosity(1);
    int i, j, isotope, nb_lev_ch3oh, nb_vibr_ch3oh, ang_mom_max;
    double t, rate, spin, mol_mass;

    string fname;
    ofstream outfile;

    fname = output_path + "ch3oha_oh2_extrap.txt";
    outfile.open(fname.c_str());
    outfile.setf(ios::scientific);
    outfile << setprecision(4);

    //  A-class methanol
    nb_lev_ch3oh = 768;
    nb_vibr_ch3oh = 2;
    ang_mom_max = 15;

    mol_mass = 32. * ATOMIC_MASS_UNIT;
    molecule ch3oh_a_mol("CH3OHa", isotope = 1, mol_mass, spin = 1.5);

    ch3oh_diagram* ch3oh_a_levels =
        new ch3oh_diagram(input_data_path, ch3oh_a_mol, nb_lev_ch3oh, nb_vibr_ch3oh, ang_mom_max, verbosity);

    ch3oh_ph2_coll_data* ch3oh_ph2_cd =
        new ch3oh_ph2_coll_data(input_data_path, ch3oh_a_levels, verbosity);

    for (i = 1; i < nb_lev_ch3oh; i++) {
        for (j = 0; j < i; j++) {
            outfile << left << setw(5) << i << setw(5) << j;

            t = 1.;
            rate = ch3oh_ph2_cd->get_rate(i, j, 100.);
            if (ch3oh_ph2_cd->get_rate(i, j, 200.) > MIN_COLLISION_RATE && ch3oh_ph2_cd->get_rate(i, j, 100.) > MIN_COLLISION_RATE)
                t = rate * sqrt(2) / ch3oh_ph2_cd->get_rate(i, j, 200.);

            outfile << t << endl;
        }
    }
    outfile.close();
}
