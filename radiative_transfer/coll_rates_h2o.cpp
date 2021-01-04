//
// Modifications:
// 19.07.2016. The class h2o_collisions was rewritten. Two data sets were added: collisions with electrons and H atom;
// 06.03.2017. Check for errors.
// 11.09.2017. Check for errors. Calculations of H2O-He rovibrational rates were not checked.

#include <stdlib.h>
#include <memory.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

#include "constants.h"
#include "utils.h"
#include "coll_rates_h2o.h"

#define MAX_TEXT_LINE_WIDTH 1000
#define SOURCE_NAME "coll_rates_h2o.cpp"
using namespace std;

//
// The classes containing collisional data
//

// The rate coefficients for the collisional transitions between lowest 45 levels of ortho(para) h2o molecule in h2o-h2 collisions; 
// A. Faure et al., Astron. Astrophys. 472, 1029 (2007).
h2o_oh2_coll_data::h2o_oh2_coll_data(const string &data_path, const energy_diagram *h2o_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, l, li, lf, nb_lines;
	string file_name;
	ifstream input;

	nb_lev = 45;
	imax = nb_lev*(nb_lev-1)/2;
	jmax = 9; // one point is reserved for 0 K;

	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));

	if (rounding(h2o_di->mol.spin) == 0) 
		file_name = data_path + "coll_h2o/coll_ph2o_oh2.txt";
	else file_name = data_path + "coll_h2o/coll_oh2o_oh2.txt";	

	input.open(file_name.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lines;
	
	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	for (i = 0; i < nb_lines && i < imax; i++)
	{
		input >> l >> li >> lf;
		for (j = 1; j < jmax; j++) {
			input >> coeff[i][j];
		}
	}
	input.close();
	calc_coeff_deriv();
	
	if (verbosity) {
		cout << left << "  data have been read from file " << file_name << endl 
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

h2o_ph2_coll_data::h2o_ph2_coll_data(const string &data_path, const energy_diagram *h2o_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, l, li, lf, nb_lines;
	string file_name;
	ifstream input;

	nb_lev = 45;
	imax = nb_lev*(nb_lev-1)/2;
	jmax = 9; // one point is reserved for 0 K;

	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	if (rounding(h2o_di->mol.spin) == 0) 
		file_name = data_path + "coll_h2o/coll_ph2o_ph2.txt";
	else file_name = data_path + "coll_h2o/coll_oh2o_ph2.txt";		

	input.open(file_name.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in :" << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lines;
	
	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	for (i = 0; i < nb_lines && i < imax; i++)
	{
		input >> l >> li >> lf;
		for (j = 1; j < jmax; j++) {
			input >> coeff[i][j];
		}
	}
	input.close();
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

// Rate coefficients for the rovibrational transitions between 411 levels of ortho-h2o and 413 levels of para-h2o in h2o-h2 collisions; 
// A. Faure and E. Josselin, Astron. Astrophys. 492, 257 (2008).
h2o_h2_coll_rovibr_data::h2o_h2_coll_rovibr_data(const string &data_path, const energy_diagram *h2o_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, v1, v2, j1, j2, tau1, tau2, up_lev, low_lev, line, nb_lines;
	double val;
	
	string file_name;
	ifstream input;

	nb_lev = h2o_di->nb_lev;	
	imax = nb_lev*(nb_lev-1)/2;
	jmax = 12; // one point is reserved for 0 K;

	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax *jmax *sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax *jmax *sizeof(double));

	if (rounding(h2o_di->mol.spin) == 0) 
		file_name = data_path + "coll_h2o/coll_ph2o_h2_rovibr.txt";
	else file_name = data_path + "coll_h2o/coll_oh2o_h2_rovibr.txt";	

	input.open(file_name.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional rovibr data " << file_name << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lines;
	
	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	for (line = 0; line < nb_lines; line++)
	{
		input >> v1 >> j1 >> tau1 >> v2 >> j2 >> tau2;

		up_lev = h2o_di->get_nb(v1, j1, tau1);
		low_lev = h2o_di->get_nb(v2, j2, tau2);

		if (low_lev != -1 && up_lev != -1)
		{
			i = up_lev*(up_lev-1)/2 + low_lev;				
			for (j = 1; j < jmax; j++) {
				input >> coeff[i][j];
			}
		}
		else {
			for (j = 1; j < jmax; j++) {
				input >> val;
			}
		}
	}
	input.close();
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

// Rate coefficients for the collisional transitions between lowest 45 levels of ortho(para) h2o in h2o-he collisions; 
// S. Green et al., Astrophys. J. Suppl. Ser. 85, 181 (1993).
h2o_he_coll_data::h2o_he_coll_data(const string &data_path, const energy_diagram *h2o_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, f, nb, li, lf, nb_lines;
	double **temp;
	
	string file_name;
	ifstream input;

	nb_lev = 45;
	imax = nb_lev*(nb_lev-1)/2;
	jmax = 11; // one point is reserved for 0 K;

	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax *jmax *sizeof(double));

	temp = alloc_2d_array<double>(2*imax, jmax);
	memset(*temp, 0, 2*imax *jmax *sizeof(double));

	if (rounding(h2o_di->mol.spin) == 0) 
		file_name = data_path + "coll_h2o/coll_ph2o_he.txt";
	else file_name = data_path + "coll_h2o/coll_oh2o_he.txt";

	input.open(file_name.c_str(), ios_base::in);
	
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lines;
	
	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	for (i = 0; i < nb_lines && i < 2*imax; i++)
	{
		input >> li >> lf >> f >> f;
		for (j = 1; j < jmax; j++) {
			input >> temp[i][j];
		}
	}
	input.close();

	for (li = 1; li < nb_lev; li++) 
	{
		for (lf = 0; lf < li; lf++)
		{
			i = li *(nb_lev-1) + lf;
			f = lf *(nb_lev-1) + li - 1;
			nb = li*(li-1)/2 + lf;

			if (li < h2o_di->nb_lev) {
				for (j = 1; j < jmax; j++)
				{
					coeff[nb][j] = 0.5*(temp[i][j] + temp[f][j] *h2o_di->lev_array[lf].g / ((double) h2o_di->lev_array[li].g)
						*exp((h2o_di->lev_array[li].energy - h2o_di->lev_array[lf].energy) *CM_INVERSE_TO_KELVINS/tgrid[j])); 
				}
			}
		}
	}
	calc_coeff_deriv();
	free_2d_array<double>(temp);

	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

// Extrapolated data of Green et al. (1993) to higher rovibrational levels of h2o, 411 levels of ortho(para) h2o;
h2o_he_coll_rovibr_data::h2o_he_coll_rovibr_data(const string &data_path, const energy_diagram *h2o_di, bool is_scaled, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, v1, v2, j1, j2, tau1, tau2, up_lev, low_lev, line, nb_lines;
	double val;
	
	string file_name;
	ifstream input;

	nb_lev = h2o_di->nb_lev;	
	imax = nb_lev*(nb_lev-1)/2;
	
	if (rounding(h2o_di->mol.spin) == 0) {
		if (is_scaled) 
			file_name = data_path + "coll_h2o/coll_ph2o_he_rovibr_scaled.txt";
		else file_name = data_path + "coll_h2o/coll_ph2o_he_rovibr.txt";
	}
	else {
		if (is_scaled) 
			file_name = data_path + "coll_h2o/coll_oh2o_he_rovibr_scaled.txt";
		else file_name = data_path + "coll_h2o/coll_oh2o_he_rovibr.txt";
	}

	input.open(file_name.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional rovibr data h2o-he " << file_name << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lines >> jmax;
	jmax++; // one point is reserved for 0 K;

	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax *jmax *sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax *jmax *sizeof(double));

	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	for (line = 0; line < nb_lines; line++)
	{
		input >> v1 >> j1 >> tau1 >> v2 >> j2 >> tau2;
		
		up_lev = h2o_di->get_nb(v1, j1, tau1);
		low_lev = h2o_di->get_nb(v2, j2, tau2);

		if (low_lev != -1 && up_lev != -1)
		{
			i = up_lev*(up_lev-1)/2 + low_lev;			
			for (j = 1; j < jmax; j++) {
				input >> coeff[i][j];
			}
		}
		else {
			for (j = 1; j < jmax; j++) {
				input >> val;
			}
		}
	}
	input.close();
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

// Rate coefficients for the rovibrational transitions between 411 levels of ortho-h2o and 413 levels of para-h2o in h2o-e collisions; 
// A. Faure and E. Josselin, Astron. Astrophys. 492, 257 (2008).
h2o_e_coll_rovibr_data::h2o_e_coll_rovibr_data(const string &data_path, const energy_diagram *h2o_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, v1, v2, j1, j2, tau1, tau2, up_lev, low_lev, line, nb_lines;
	double val;
	
	string file_name;
	ifstream input;

	nb_lev = h2o_di->nb_lev;	
	imax = nb_lev*(nb_lev-1)/2;
	jmax = 12; // one point is reserved for 0 K;

	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax *jmax *sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax *jmax *sizeof(double));

	if (rounding(h2o_di->mol.spin) == 0) 
		file_name = data_path + "coll_h2o/coll_ph2o_e_rovibr.txt";
	else file_name = data_path + "coll_h2o/coll_oh2o_e_rovibr.txt";	

	input.open(file_name.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional rovibr data " << file_name << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lines;
	
	tgrid[0] = 0.; // point at 0 K;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	for (line = 0; line < nb_lines; line++)
	{
		input >> v1 >> j1 >> tau1 >> v2 >> j2 >> tau2;

		up_lev = h2o_di->get_nb(v1, j1, tau1);
		low_lev = h2o_di->get_nb(v2, j2, tau2);

		if (low_lev != -1 && up_lev != -1)
		{
			i = up_lev*(up_lev-1)/2 + low_lev;				
			for (j = 1; j < jmax; j++) {
				input >> coeff[i][j];
			}
		}
		else {
			for (j = 1; j < jmax; j++) {
				input >> val;
			}
		}
	}
	input.close();
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

// Rate coefficients for transitions between 45 levels of ortho/para-H2O and H atom
// Daniel et al., MNRAS 446, p. 2312, 2015;
h2o_h_coll_data::h2o_h_coll_data(const string &path, const energy_diagram *h2o_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, l, li, lf, nb_lines;
	string file_name;
	ifstream input;

	nb_lev = 45;
	imax = nb_lev*(nb_lev-1)/2;
	jmax = 15; // one point is reserved for 0 K;

	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	if (rounding(h2o_di->mol.spin) == 0) 
		file_name = path + "coll_h2o/coll_ph2o_h.txt";
	else file_name = path + "coll_h2o/coll_oh2o_h.txt";	

	input.open(file_name.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lines;
	
	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	for (i = 0; i < nb_lines && i < imax; i++)
	{
		input >> li >> lf >> l >> l;
		for (j = 1; j < jmax; j++) {
			input >> coeff[i][j];
		}
	}
	input.close();
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

//
// The class calculates collisional rates
//

h2o_collisions::h2o_collisions(const string &data_path, const energy_diagram* h2o_di, bool he_is_scaled, int verbosity)
{
	if (verbosity) 
		cout << "H2O collisional rate coefficients are being initializing..." << endl;
	
	nb_lev = h2o_di->nb_lev;

	// data for rovibrational transitions includes temperatures >= 200 K (for H2) and >= 120 K (for He);
	coll_data.push_back( new h2o_he_coll_data(data_path, h2o_di, verbosity) );
	coll_data.push_back( new h2o_he_coll_rovibr_data(data_path, h2o_di, he_is_scaled, verbosity) );
	coll_data.push_back( new h2o_ph2_coll_data(data_path, h2o_di, verbosity) );
	coll_data.push_back( new h2o_oh2_coll_data(data_path, h2o_di, verbosity) );
	coll_data.push_back( new h2o_h2_coll_rovibr_data(data_path, h2o_di, verbosity) );
	coll_data.push_back( new h2o_h_coll_data(data_path, h2o_di, verbosity) );
	nb1 = (int) coll_data.size();

	// data on electron collisions must be here;
	coll_data.push_back( new h2o_e_coll_rovibr_data(data_path, h2o_di, verbosity) );
	nb2 = (int) coll_data.size();

	// data on H+ collisions must be here;
	nb3 = (int) coll_data.size();

	max_temp = new double [nb3];
	for (int i = 0; i < nb3; i++) {
		max_temp[i] = coll_data[i]->get_max_temp();
	}
}

void h2o_collisions::set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, 
	double h_conc, double el_conc, double *&concentration, int *&indices) const
{
	collisional_transitions::set_gas_param(temp_neutrals, temp_el, he_conc, ph2_conc, oh2_conc, h_conc, el_conc, concentration, 
		indices);

	concentration[0] = he_conc;
	concentration[1] = he_conc + 0.2*h_conc; // collisions with H atoms are taken into account;
	concentration[2] = ph2_conc;
	concentration[3] = oh2_conc;
	concentration[4] = ph2_conc + oh2_conc;
	concentration[5] = h_conc;
}

// The energy of the first level is higher, up_lev.nb > low_lev.nb;
void h2o_collisions::get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, 
	double &up_rate, double temp_neutrals, const double *concentration, const int *indices) const
{
	if (up_lev.nb < 45) {
		down_rate = coll_data[0]->get_rate(up_lev.nb, low_lev.nb, indices[0], (temp_neutrals < max_temp[0]) ? temp_neutrals : max_temp[0]) *concentration[0]
			+ coll_data[2]->get_rate(up_lev.nb, low_lev.nb, indices[2], (temp_neutrals < max_temp[2]) ? temp_neutrals : max_temp[2]) *concentration[2] 
			+ coll_data[3]->get_rate(up_lev.nb, low_lev.nb, indices[3], (temp_neutrals < max_temp[3]) ? temp_neutrals : max_temp[3]) *concentration[3]
			+ coll_data[5]->get_rate(up_lev.nb, low_lev.nb, indices[5], (temp_neutrals < max_temp[5]) ? temp_neutrals : max_temp[5]) *concentration[5];
	}
	else {
		// there is no upper restriction on the level nb for these data;
		down_rate = coll_data[1]->get_rate(up_lev.nb, low_lev.nb, indices[1], (temp_neutrals < max_temp[1]) ? temp_neutrals : max_temp[1]) *concentration[1] 
			+ coll_data[4]->get_rate(up_lev.nb, low_lev.nb, indices[4], (temp_neutrals < max_temp[4]) ? temp_neutrals : max_temp[4]) *concentration[4];
	}

	if (down_rate > MIN_COLLISION_RATE)
		up_rate = down_rate *exp((low_lev.energy - up_lev.energy)*CM_INVERSE_TO_KELVINS/temp_neutrals) *up_lev.g / ((double) low_lev.g);
	else up_rate = down_rate = 0.;
}

//
// Functions
//

// Check if files with original data are present; 
void reformat_h2o_coll_data(const std::string &path, int spin)
{
	int i, j, k, jmax, up_lev, low_lev, line, nb_lines, nb_lev;
	double val, mass;

	string fname_in[2], fname_out[2], pname[2];
	ifstream input;
	ofstream output;

	pname[0] = "electrons";
	pname[1] = "H2";

	if (spin == 0) 
	{
		nb_lev = 413;
		fname_in[0] = path + "coll_h2o/ph2o_e_faure2008.txt";
		fname_out[0] = path + "coll_h2o/coll_ph2o_e_rovibr.txt";

		fname_in[1] = path + "coll_h2o/coll-ph2o-h2_faure2008.txt";
		fname_out[1] = path + "coll_h2o/coll_ph2o_h2_rovibr.txt";
	}
	else 
	{
		spin = 1;
		nb_lev = 411;
		fname_in[0] = path + "coll_h2o/oh2o_e_faure2008.txt";
		fname_out[0] = path + "coll_h2o/coll_oh2o_e_rovibr.txt";

		fname_in[1] = path + "coll_h2o/coll-oh2o-h2_faure2008.txt";
		fname_out[1] = path + "coll_h2o/coll_oh2o_h2_rovibr.txt";
	}
	
	mass = 18.*ATOMIC_MASS_UNIT;
	molecule h2o_mol("h2o", 1, mass, spin);

	h2o_diagram *h2o_di =
		new h2o_diagram(path, h2o_mol, nb_lev);

	for (k = 0; k < 2; k++)
	{
		output.open(fname_out[k].c_str(), ios_base::out);
		input.open(fname_in[k].c_str(), ios_base::in);

		if (!input) {
			cout << "Error in " << SOURCE_NAME << ": can't open file with collisional rovibr data " << fname_in << endl;
			exit(1);
		}

		output << "# Faure A., Josselin E., Astron. Astrophys. 492, pp. 257-264, 2008, BASECOL V1 P2012-06" << endl;
		if (spin == 0)
			output << "# collisions of para-H2O with " << pname[k] << ", 413 levels" << endl;
		else output << "# collisions of ortho-H2O with " << pname[k] << ", 411 levels" << endl;
	
		// there is no comments in the input file,
		input >> nb_lines;
		output << nb_lines << endl;
	
		jmax = 11; 
		for (j = 0; j < jmax; j++) 
		{
			input >> val;
			output << left << setw(12) << val;
		}
		output << endl;

		for (line = 0; line < nb_lines; line++)
		{
			if (k == 0)
				input >> up_lev >> low_lev >> i >> i;
			else input >> i >> up_lev >> low_lev;

			up_lev--;
			low_lev--;

			output << left << setw(5) << h2o_di->lev_array[up_lev].v << setw(5) << rounding(h2o_di->lev_array[up_lev].j) 
				<< setw(5) << rounding(h2o_di->lev_array[up_lev].k1 - h2o_di->lev_array[up_lev].k2) 
				<< setw(5) << h2o_di->lev_array[low_lev].v << setw(5) << rounding(h2o_di->lev_array[low_lev].j) 
				<< setw(5) << rounding(h2o_di->lev_array[low_lev].k1 - h2o_di->lev_array[low_lev].k2);

			for (j = 0; j < jmax; j++) 
			{
				input >> val;
				output << left << setw(12) << val;
			}
			if (line < nb_lines-1)
				output << endl;
		}
		input.close();
		output.close();
	}
}
