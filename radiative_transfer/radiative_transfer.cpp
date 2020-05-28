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
#include "cloud_data.h"
#include "iteration_lvg.h"
#include "iteration_control.h"
#include "transition_data.h"
#include "transition_list.h"
#include "maser_luminosity.h"

#define NB_OF_OMP_PROCESSES 2 // not used
#define NB_JOINED_LAYERS 2
#define MAX_NB_ITER_OH 20000

// These parameters must be the same as in the C-type shock simulations (are necessary for dust model),
#define HE_TO_H_NB_RATIO 0.09
#define STANDARD_NB_CR_PHOTONS 3000.
#define MICROTURBULENT_SPEED 3.e+4 // in cm/s, must be the same as in the shock simulations,

using namespace std;

// the relative population error, is used in calculations of populations and in finding of inverted transitions,
const double rel_population_error = 1.e-5;
//const double population_error = 1.e-8; // !!!

// OH radiative transfer calculations - with line overlap,
void calc_oh_masers(string input_data_path, string sim_data_path, string output_path, bool save_line_stat);

// OH radiative transfer calculations - no line overlap,
void calc_oh_masers_test(string input_data_path, string sim_data_path, string output_path);

// CH3OH radiative transfer calculations - with line overlap,
void calc_ch3oh_masers(string input_data_path, string sim_data_path, string output_path, bool save_line_stat);

// the calculation of the line overlap number as a function of doppler width,
void calc_nb_line_overlaps(string input_data_path, string output_path);

void calc_molecular_populations(cloud_data* cloud, iteration_scheme_lvg* it_scheme_lvg, 
    energy_diagram* mol_levels, einstein_coeff* mol_einst, collisional_transitions* mol_coll, double* mol_popul, int nb_lev, 
    int iter_max_nb, bool acceleration, int verbosity);

// saving the level populations in the format of the c-type shock simulations
void save_mol_data(std::string fname, cloud_data* cloud, double *lev_pop, int nb_lev);

int main()
{
    bool save_line_stat;
    int i;
    string sim_data_path, input_data_path, output_path;
    stringstream sstr;

    // the path to the directory with spectroscopic data,
    input_data_path = "C:/Users/Александр/Documents/input_data/";
    
//    calc_nb_line_overlaps(input_data_path, output_path = "");

    //sim_data_path = "C:/Users/Александр/Documents/Данные и графики/paper C-type shocks - new data on H-H2 collisions/output_data_2e4_magnf-b2/shock_cr1-15_25_new/";
    sim_data_path = "C:/Users/Александр/Documents/Данные и графики/paper C-type shocks - new data on H-H2 collisions/output_data_2e5/shock_cr1-15_20/";
    
    output_path = sim_data_path + "radiative_transfer/masers/";
//    calc_oh_masers(input_data_path, sim_data_path, output_path, save_line_stat = true); // with line statistics
    
    output_path = sim_data_path + "radiative_transfer/masers_oh_test/";
//    calc_oh_masers_test(input_data_path, sim_data_path, output_path);

    sim_data_path = "C:/Users/Александр/Documents/Данные и графики/paper C-type shocks - new data on H-H2 collisions/output_data_2e3_magnf-b2/shock_cr1-15_";
//    sim_data_path = "C:/Users/Александр/Documents/Данные и графики/paper C-type shocks - new data on H-H2 collisions/output_data_2e4/shock_cr1-15_";
    
    for (i = 5; i <= 100; i += 5) {
        sstr.clear();
        sstr.str("");

        if (i < 10) 
            sstr << sim_data_path << "0" << i << "/";
        else 
            sstr << sim_data_path << i << "/";
        
        cout << sstr.str() << endl;
        output_path = sstr.str() + "radiative_transfer/masers/";
//        output_path = sstr.str() + "radiative_transfer/masers_oh_coll_ext/";
        calc_oh_masers(input_data_path, sstr.str(), output_path, save_line_stat = false); // no line statistics

//        output_path = sstr.str() + "radiative_transfer/masers_oh_test/";
//        calc_oh_masers_test(input_data_path, sstr.str(), output_path);
    }
}

void calc_oh_masers(string input_data_path, string sim_data_path, string output_path, bool save_line_stat)
{
    bool acceleration, do_simdata_exist;
	int verbosity = 1;
	int i, isotope, nb_lev_oh, nb_of_dust_comp, nb_cloud_lay, iter_max_nb;
	
	double mol_mass, spin, dg_ratio, c_abund_pah;
	double * oh_popul;

    stringstream sstr;
    list<transition> trans_list;
    list<transition_data>::const_iterator td_it;

#ifdef _OPENMP
	omp_set_num_threads(NB_OF_OMP_PROCESSES);

#pragma omp parallel 
	{
#pragma omp master 
		{
			cout << "OpenMP is supported" << endl;
			cout << "Nb of threads: " << omp_get_num_threads() << endl;
		}
	}
#endif

	// dust model:
	c_abund_pah = 0.;
	dust_model *dust = 
		new two_component_dust_model(input_data_path, c_abund_pah, dg_ratio = 0.01, HE_TO_H_NB_RATIO, STANDARD_NB_CR_PHOTONS, verbosity);
	
	nb_of_dust_comp = dust->nb_of_comp;

	cloud_layer clayer;
	cloud_data* cloud
		= new cloud_data();
	
    // physical parameters of the particular shock model 
    do_simdata_exist = set_physical_parameters(sim_data_path, cloud);
    if (!do_simdata_exist)
        return;

    do_simdata_exist = set_molecular_conc(sim_data_path, "OH", cloud);
    if (!do_simdata_exist)
        return;

	join_layers(cloud, NB_JOINED_LAYERS);
	nb_cloud_lay = cloud->nb_lay;

	cloud->set_vel_turb(MICROTURBULENT_SPEED);
	cloud->set_dust_model(dust);
	
    // OH with hyperfine levels, 24, 56
    // up to the HF 20 levels may be involved in the pumping of the maser (Gray, Proc. IAU Symp. 287, 2012)  
    nb_lev_oh = 56;
    oh_popul = new double[nb_cloud_lay * nb_lev_oh];

    mol_mass = 17. * ATOMIC_MASS_UNIT;
    molecule oh_mol("OH", isotope = 1, mol_mass, spin = 0.5);

    oh_hf_diagram *oh_levels = 
        new oh_hf_diagram(input_data_path, oh_mol, nb_lev_oh, verbosity);
    
    oh_hf_einstein_coeff *oh_einst = 
        new oh_hf_einstein_coeff(input_data_path, oh_levels, verbosity);

    oh_hf_collisions *oh_coll = 
        new oh_hf_collisions(input_data_path, oh_levels, verbosity);

    iteration_scheme_line_overlap * it_scheme_loverlap 
        = new iteration_scheme_line_overlap(input_data_path, verbosity);

    it_scheme_loverlap->init_molecule_data(oh_levels, oh_einst, oh_coll);
    it_scheme_loverlap->set_dust_model(dust);

    calc_molecular_populations(cloud, it_scheme_loverlap, oh_levels, oh_einst, oh_coll, oh_popul, nb_lev_oh,
        iter_max_nb = MAX_NB_ITER_OH, acceleration = false, verbosity = false);

    transition_data_container* trans_data_oh
        = new transition_data_container(cloud, oh_levels, oh_einst);

    trans_data_oh->find(oh_popul, rel_population_error); 
    
    // adding transitions,
    trans_list.clear();
    oh_trans_list(oh_levels, trans_list);
    trans_data_oh->add(trans_list, oh_popul);

    lim_luminosity_lvg(it_scheme_loverlap, trans_data_oh, cloud, oh_levels, oh_einst, oh_coll, oh_popul);
    
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
    sstr << output_path + "radtr_data_" + oh_mol.name + ".txt";
    save_mol_data(sstr.str(), cloud, oh_popul, nb_lev_oh);

    // line statistic for cloud layers in the post shock gas:
    if (save_line_stat) {
        for (i = nb_cloud_lay-1; i > nb_cloud_lay/2; i -= 5) {
            clayer = cloud->lay_array[i];
            it_scheme_loverlap->set_vel_grad(clayer.velg_n);
            it_scheme_loverlap->set_dust_parameters(clayer.dust_grain_conc, clayer.dust_grain_temp);
            it_scheme_loverlap->set_parameters(clayer.temp_n, clayer.temp_el, clayer.el_conc, clayer.h_conc, clayer.ph2_conc, clayer.oh2_conc, clayer.he_conc,
                clayer.mol_conc, clayer.vel_turb);

            sstr.clear();
            sstr.str("");
            sstr << output_path + oh_mol.name + "_line_stat_" << i << ".txt";
            it_scheme_loverlap->calc_line_stat(sstr.str(), oh_popul+ i * nb_lev_oh);
        }
    }

    delete[] oh_popul;
    delete oh_levels;
    delete oh_einst;
    delete oh_coll;
    delete trans_data_oh;
    delete it_scheme_loverlap;
    delete dust;
    delete cloud;
}

void calc_oh_masers_test(string input_data_path, string sim_data_path, string output_path)
{
    bool acceleration, do_simdata_exist;
    int verbosity = 1;
    int nb_of_dust_comp, nb_cloud_lay, nb_lev_oh, isotope, iter_max_nb;
    
    double c_abund_pah, dg_ratio, mol_mass, spin;
    double* oh_popul;

    stringstream sstr;
    list<transition> trans_list;
    list<transition_data>::const_iterator td_it;

    // dust model:
    c_abund_pah = 0.;
    dust_model* dust =
        new two_component_dust_model(input_data_path, c_abund_pah, dg_ratio = 0.01, HE_TO_H_NB_RATIO, STANDARD_NB_CR_PHOTONS, verbosity);

    nb_of_dust_comp = dust->nb_of_comp;

    // physical parameters of the particular shock model 
    cloud_data* cloud
        = new cloud_data();

    do_simdata_exist = set_physical_parameters(sim_data_path, cloud);
    if (!do_simdata_exist)
        return;

    do_simdata_exist = set_molecular_conc(sim_data_path, "OH", cloud);
    if (!do_simdata_exist)
        return;

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

    iteration_scheme_lvg* it_scheme_lvg
        = new iteration_scheme_lvg(input_data_path, verbosity);

    it_scheme_lvg->set_dust_model(dust);
    it_scheme_lvg->init_molecule_data(oh_levels, oh_einst, oh_coll);

    calc_molecular_populations(cloud, it_scheme_lvg, oh_levels, oh_einst, oh_coll, oh_popul, nb_lev_oh, 
        iter_max_nb = MAX_NB_ITER_OH, acceleration = false, verbosity = false);

    transition_data_container* trans_data_oh
        = new transition_data_container(cloud, oh_levels, oh_einst);
 
    trans_data_oh->find(oh_popul, rel_population_error);
       
    trans_list.clear();
    oh_trans_list(oh_levels, trans_list);
    trans_data_oh->add(trans_list, oh_popul);

    lim_luminosity_lvg(it_scheme_lvg, trans_data_oh, cloud, oh_levels, oh_einst, oh_coll, oh_popul);
    
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
    delete cloud;
}

void calc_molecular_populations(cloud_data* cloud, iteration_scheme_lvg* it_scheme_lvg, 
    energy_diagram *mol_levels, einstein_coeff *mol_einst, collisional_transitions *mol_coll, double *mol_popul, int nb_lev, 
    int iter_max_nb, bool acceleration, int verbosity)
{
    bool is_solution_found, is_solution_found_prev;
    int i, lay_nb, nb_cloud_lay;
    
    vector<int> bad_layers;
    cloud_layer clayer;

    iteration_control<iteration_scheme_lvg> it_control(it_scheme_lvg);

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
            is_solution_found =
                it_control.calculate_populations(mol_popul + lay_nb * nb_lev, iter_max_nb, rel_population_error, acceleration, verbosity);
        }

        if (!is_solution_found) {
            boundary_layer_populations(mol_popul + lay_nb * nb_lev, mol_levels, mol_einst, mol_coll,
                clayer.temp_n, clayer.temp_el, clayer.el_conc, clayer.h_conc, clayer.ph2_conc, clayer.oh2_conc, clayer.he_conc);

            is_solution_found =
                it_control.calculate_populations(mol_popul + lay_nb * nb_lev, iter_max_nb, rel_population_error, acceleration, verbosity);
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

void save_mol_data(std::string fname, cloud_data* cloud, double* lev_pop, int nb_lev)
{
	int i, j;
	double a;
	ofstream output;

	output.open(fname.c_str());
	output << scientific;
	output.precision(4);
	output << left << "! three parameters are given: shock speed (cm/s), turbulent velocity (cm/s), number of specimen levels:" << endl
		<< setw(13) << cloud->lay_array[0].vel_n << setw(13) << cloud->lay_array[0].vel_turb << nb_lev << endl;

	output << left << setw(18) << "!depth(cm)" << setw(13) << "gas_temp(K)" << setw(13) << "el_temp" << setw(13) << "dust_temp(K)"
		<< setw(13) << "gasvel(cm/s)" << setw(13) << "concHe(cm-3)" << setw(13) << "conc_pH2" << setw(13) << "conc_oH2"
		<< setw(13) << "conc_H" << setw(13) << "conc_e" << setw(13) << "conc_H_nucl" << setw(13) << "conc_mol";
	
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
			a = lev_pop[i*nb_lev + j];
			if (a < 1.e-99)
				a = 0.;
			output << left << setw(13) << a* cloud->lay_array[i].mol_conc;
		}
	}
}

void calc_ch3oh_masers(string input_data_path, string sim_data_path, string output_path, bool save_line_stat)
{
}

/*
bool acceleration, do_simdata_exist;
    int verbosity = 1;
    int i, isotope, nb_lev_ch3oh, nb_vibr_ch3oh, ang_mom_max, nb_lev_onh3, nb_lev_pnh3, nb_lev_oh, nb_of_dust_comp, nb_cloud_lay, iter_max_nb;

    double f, vel_turb, mol_mass, spin, dg_ratio, c_abund_pah, temp_stat;
    double *ch3oh_a_popul, * ch3oh_e_popul, * onh3_popul, * pnh3_popul, * oh_popul;

       // no line overlap,
    iteration_scheme_lvg* it_scheme_lvg
        = new iteration_scheme_lvg(input_data_path, verbosity);

    it_scheme_lvg->set_dust_model(dust);

    iteration_control<iteration_scheme_lvg> it_control(it_scheme_lvg);
    it_control.set_accel_parameters(25, 4, 5);

#if (!CALC_ONLY_OH)
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

do_simdata_exist = set_molecular_conc(sim_data_path, "CH3OH", cloud, f = 0.5);
if (!do_simdata_exist)
return;

it_scheme_lvg->init_molecule_data(ch3oh_a_levels, ch3oh_a_einst, ch3oh_a_coll);
calc_molecular_populations(cloud, it_scheme_lvg, it_control, ch3oh_a_levels, ch3oh_a_einst, ch3oh_a_coll, ch3oh_a_popul, nb_lev_ch3oh);

transition_data_container* trans_data_ch3oh_a
= new transition_data_container(cloud, ch3oh_a_levels, ch3oh_a_einst);

// only transitions with inversion of level populations are saved individually,
trans_data_ch3oh_a->find(ch3oh_a_popul, rel_population_error);
lim_luminosity_lvg(it_scheme_lvg, trans_data_ch3oh_a, cloud, ch3oh_a_levels, ch3oh_a_einst, ch3oh_a_coll, ch3oh_a_popul);

trans_list.clear();
ch3oh_classI_trans_list(ch3oh_a_levels, trans_list);
//ch3oh_classII_trans_list(ch3oh_a_levels, trans_list);

// only those transitions are saved, that are in the given list
td_it = trans_data_ch3oh_a->data.begin();
while (td_it != trans_data_ch3oh_a->data.end())
{
    tl_it = find(trans_list.begin(), trans_list.end(), *(td_it->trans));
    if (tl_it != trans_list.end()) {
        sstr.clear();
        sstr.str("");
        sstr << output_path + "inv_trans_" + ch3oh_a_mol.name + "_";
        sstr << (int)(td_it->trans->freq * 1.e-6);
        sstr << ".txt";
        trans_data_ch3oh_a->save_transition(sstr.str(), *(td_it->trans));
    }
    td_it++;
}
// adding more transitions:
// the limiting luminosity is not calculated for these transitions,
trans_data_ch3oh_a->add(trans_list, ch3oh_a_popul);

sstr.clear();
sstr.str("");
sstr << output_path + "inv_trans_" << ch3oh_a_mol.name << ".txt";
trans_data_ch3oh_a->save_full_data(sstr.str());

sstr.clear();
sstr.str("");
sstr << output_path + "radtr_data_" << ch3oh_a_mol.name << ".txt";
save_mol_data(sstr.str(), cloud, ch3oh_a_popul, nb_lev_ch3oh);

delete[] ch3oh_a_popul;
delete ch3oh_a_levels;
delete ch3oh_a_einst;
delete ch3oh_a_coll;
delete trans_data_ch3oh_a;

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

do_simdata_exist = set_molecular_conc(sim_data_path, "CH3OH", cloud, f = 0.5);
if (!do_simdata_exist)
return;

it_scheme_lvg->init_molecule_data(ch3oh_e_levels, ch3oh_e_einst, ch3oh_e_coll);
calc_molecular_populations(cloud, it_scheme_lvg, it_control, ch3oh_e_levels, ch3oh_e_einst, ch3oh_e_coll, ch3oh_e_popul, nb_lev_ch3oh);

transition_data_container* trans_data_ch3oh_e
= new transition_data_container(cloud, ch3oh_e_levels, ch3oh_e_einst);

trans_data_ch3oh_e->find(ch3oh_e_popul, rel_population_error);
lim_luminosity_lvg(it_scheme_lvg, trans_data_ch3oh_e, cloud, ch3oh_e_levels, ch3oh_e_einst, ch3oh_e_coll, ch3oh_e_popul);

trans_list.clear();
ch3oh_classI_trans_list(ch3oh_e_levels, trans_list);
//ch3oh_classII_trans_list(ch3oh_e_levels, trans_list);

td_it = trans_data_ch3oh_e->data.begin();
while (td_it != trans_data_ch3oh_e->data.end())
{
    tl_it = find(trans_list.begin(), trans_list.end(), *(td_it->trans));
    if (tl_it != trans_list.end()) {
        sstr.clear();
        sstr.str("");
        sstr << output_path + "inv_trans_" + ch3oh_e_mol.name + "_";
        sstr << (int)(td_it->trans->freq * 1.e-6);
        sstr << ".txt";
        trans_data_ch3oh_e->save_transition(sstr.str(), *(td_it->trans));
    }
    td_it++;
}
trans_data_ch3oh_e->add(trans_list, ch3oh_e_popul);

sstr.clear();
sstr.str("");
sstr << output_path + "inv_trans_" << ch3oh_e_mol.name << ".txt";
trans_data_ch3oh_e->save_full_data(sstr.str());

sstr.clear();
sstr.str("");
sstr << output_path + "radtr_data_" << ch3oh_e_mol.name << ".txt";
save_mol_data(sstr.str(), cloud, ch3oh_e_popul, nb_lev_ch3oh);

delete[] ch3oh_e_popul;
delete ch3oh_e_levels;
delete ch3oh_e_einst;
delete ch3oh_e_coll;
delete trans_data_ch3oh_e;

// ortho-NH3
// NH3 molecule data, o-NH3 has k = 3n, n is an integer, for p-NH3 k != 3n
nb_lev_onh3 = 17; // ortho-NH3: He coll data - 22, H2 coll data - 17
onh3_popul = new double[nb_cloud_lay * nb_lev_onh3];

mol_mass = 17. * ATOMIC_MASS_UNIT;
molecule onh3_mol("oNH3", isotope = 1, mol_mass, spin = 1.5);

nh3_diagram* onh3_levels =
new nh3_diagram(input_data_path, onh3_mol, nb_lev_onh3, verbosity);

nh3_einstein_coeff* onh3_einst =
new nh3_einstein_coeff(input_data_path, onh3_levels, verbosity);

nh3_collisions* onh3_coll =
new nh3_collisions(input_data_path, onh3_levels, verbosity);

do_simdata_exist = set_molecular_conc(sim_data_path, "NH3", cloud, f = 0.5);
if (!do_simdata_exist)
return;

it_scheme_lvg->init_molecule_data(onh3_levels, onh3_einst, onh3_coll);
calc_molecular_populations(cloud, it_scheme_lvg, it_control, onh3_levels, onh3_einst, onh3_coll, onh3_popul, nb_lev_onh3);

transition_data_container* trans_data_onh3
= new transition_data_container(cloud, onh3_levels, onh3_einst);

trans_data_onh3->find(onh3_popul, rel_population_error);
lim_luminosity_lvg(it_scheme_lvg, trans_data_onh3, cloud, onh3_levels, onh3_einst, onh3_coll, onh3_popul);

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

trans_list.clear();
nh3_trans_list(onh3_levels, trans_list);
trans_data_onh3->add(trans_list, onh3_popul);

sstr.clear();
sstr.str("");
sstr << output_path + "inv_trans_" << onh3_mol.name << ".txt";
trans_data_onh3->save_full_data(sstr.str());

sstr.clear();
sstr.str("");
sstr << output_path + "radtr_data_" << onh3_mol.name << ".txt";
save_mol_data(sstr.str(), cloud, onh3_popul, nb_lev_onh3);

delete[] onh3_popul;
delete onh3_levels;
delete onh3_einst;
delete onh3_coll;
delete trans_data_onh3;

// para-NH3
nb_lev_pnh3 = 34; // para-NH3: He coll data - 16, H2 coll data - 34
pnh3_popul = new double[nb_cloud_lay * nb_lev_pnh3];

molecule pnh3_mol("pNH3", isotope = 1, mol_mass, spin = 0.5);

nh3_diagram* pnh3_levels =
new nh3_diagram(input_data_path, pnh3_mol, nb_lev_pnh3, verbosity);

nh3_einstein_coeff* pnh3_einst =
new nh3_einstein_coeff(input_data_path, pnh3_levels, verbosity);

nh3_collisions* pnh3_coll =
new nh3_collisions(input_data_path, pnh3_levels, verbosity);

do_simdata_exist = set_molecular_conc(sim_data_path, "NH3", cloud, f = 0.5);
if (!do_simdata_exist)
return;

it_scheme_lvg->init_molecule_data(pnh3_levels, pnh3_einst, pnh3_coll);
calc_molecular_populations(cloud, it_scheme_lvg, it_control, pnh3_levels, pnh3_einst, pnh3_coll, pnh3_popul, nb_lev_pnh3);

transition_data_container* trans_data_pnh3
= new transition_data_container(cloud, pnh3_levels, pnh3_einst);

// the limiting luminosity is calculated only for inverted transitions,
trans_data_pnh3->find(pnh3_popul, rel_population_error);
lim_luminosity_lvg(it_scheme_lvg, trans_data_pnh3, cloud, pnh3_levels, pnh3_einst, pnh3_coll, pnh3_popul);

trans_list.clear();
nh3_trans_list(pnh3_levels, trans_list);
trans_data_pnh3->add(trans_list, pnh3_popul);

sstr.clear();
sstr.str("");
sstr << output_path + "inv_trans_" << pnh3_mol.name << ".txt";
trans_data_pnh3->save_full_data(sstr.str());

sstr.clear();
sstr.str("");
sstr << output_path + "radtr_data_" + pnh3_mol.name + ".txt";
save_mol_data(sstr.str(), cloud, pnh3_popul, nb_lev_pnh3);

delete[] pnh3_popul;
delete pnh3_levels;
delete pnh3_einst;
delete pnh3_coll;
delete trans_data_pnh3;
#endif
*/
/*
void calc_oh_no_hfs(string input_data_path, string sim_data_path, string output_path);
void calc_oh_no_hfs(string input_data_path, string sim_data_path, string output_path)
{
    int verbosity = 1;
    int nb_of_dust_comp, nb_cloud_lay, nb_lev_oh, isotope;

    double c_abund_pah, dg_ratio, vel_turb, mol_mass, spin;
    double* oh_popul;

    stringstream sstr;
    list<transition> trans_list;
    list<transition>::const_iterator it;

    // dust model:
    c_abund_pah = 0.;
    dust_model* dust =
        new two_component_dust_model(input_data_path, c_abund_pah, dg_ratio = 0.01, HE_TO_H_NB_RATIO, STANDARD_NB_CR_PHOTONS, verbosity);

    nb_of_dust_comp = dust->nb_of_comp;

    // physical parameters of the particular shock model 
    cloud_layer clayer;
    cloud_data* cloud
        = new cloud_data();

    do_simdata_exist = set_physical_parameters(sim_data_path, cloud);
   
    join_layers(cloud, 5);
    nb_cloud_lay = cloud->nb_lay;

    // in cm/s,
    vel_turb = 5.e+4;
    cloud->set_vel_turb(vel_turb);
    cloud->set_dust_model(dust);

    iteration_scheme_lvg* it_scheme_lvg
        = new iteration_scheme_lvg(input_data_path, verbosity);

    it_scheme_lvg->set_dust_model(dust);

    iteration_control<iteration_scheme_lvg> it_control(it_scheme_lvg);
    it_control.set_accel_parameters(40, 5, 5);

    // OH with no hyperfine levels, 20
    nb_lev_oh = 20;
    oh_popul = new double[nb_cloud_lay * nb_lev_oh];

    mol_mass = 17. * ATOMIC_MASS_UNIT;
    molecule oh_mol("OH", isotope = 1, mol_mass, spin = 0.5);

    oh_diagram* oh_levels =
        new oh_diagram(input_data_path, oh_mol, nb_lev_oh, verbosity);

    oh_einstein_coeff* oh_einst =
        new oh_einstein_coeff(input_data_path, oh_levels, verbosity);

    oh_collisions* oh_coll =
        new oh_collisions(input_data_path, oh_levels, verbosity);

    set_molecular_conc(sim_data_path, "OH", cloud);
    it_scheme_lvg->init_molecule_data(oh_levels, oh_einst, oh_coll);

    calc_molecular_populations(cloud, it_scheme_lvg, it_control, oh_levels, oh_einst, oh_coll, oh_popul, nb_lev_oh);

    transition_data_container* trans_data_oh
        = new transition_data_container(cloud, oh_levels, oh_einst);

    trans_data_oh->find(oh_popul, rel_population_error);

    trans_list.clear();
    oh_trans_list(oh_levels, trans_list);
    trans_data_oh->add(trans_list, oh_popul);

    it = trans_list.begin();
    while (it != trans_list.end())
    {
        sstr.clear();
        sstr.str("");
        sstr << output_path + "inv_trans_" + oh_mol.name + "_";
        sstr << (int)(it->freq * 1.e-6);
        sstr << ".txt";
        trans_data_oh->save_transition(sstr.str(), (*it));
        it++;
    }

    sstr.clear();
    sstr.str("");
    sstr << output_path + "inv_trans_" << oh_mol.name << ".txt";
    trans_data_oh->save_full_data(sstr.str());

    sstr.clear();
    sstr.str("");
    sstr << output_path + "radtr_data_" + oh_mol.name + ".txt";
    save_mol_data(sstr.str(), cloud, oh_popul, nb_lev_oh);

    delete [] oh_popul;

    delete oh_levels;
    delete oh_einst;
    delete oh_coll;
    delete trans_data_oh;

    delete it_scheme_lvg;
    delete dust;
    delete cloud;
}*/
