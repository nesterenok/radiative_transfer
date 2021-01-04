
// 10.05.2018
//
//

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory.h>
#include <cmath>
#include <cfloat>

#include "coll_rates_nh3.h"
#include "utils.h"
#include "constants.h"

#define MAX_TEXT_LINE_WIDTH 240
#define SOURCE_NAME "coll_rates_nh3.cpp"
using namespace std;

nh3_he_coll_data::nh3_he_coll_data(const string path, const energy_diagram *di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j;
	string fname;
	ifstream input;

	if (rounding(2.*di->mol.spin) == 3) 
		fname = path + "coll_nh3/coll_onh3_he.txt";
	else fname = path + "coll_nh3/coll_pnh3_he.txt";

	input.open(fname.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lev >> jmax;
	jmax++;      // one point is reserved for 0 K;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));
	
	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	for (i = 0; i < imax; i++)
	{
		input >> j >> j >> j >> j;
		for (j = 1; j < jmax; j++) {
			input >> coeff[i][j];
		}
	}
	input.close();
	calc_coeff_deriv();
	
	if (verbosity) {
		cout << "  data have been read from file " << fname << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

nh3_ph2_coll_data::nh3_ph2_coll_data(const string path, const energy_diagram *di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j;
	string fname;
	ifstream input;

	if (rounding(2.*di->mol.spin) == 3) 
		fname = path + "coll_nh3/coll_onh3_ph2.txt";
	else fname = path + "coll_nh3/coll_pnh3_ph2.txt";

	input.open(fname.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lev >> jmax;
	jmax++;      // one point is reserved for 0 K;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));
	
	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	for (i = 0; i < imax; i++)
	{
		input >> j >> j >> j;
		for (j = 1; j < jmax; j++) {
			input >> coeff[i][j];
		}
	}
	input.close();
	calc_coeff_deriv();
	
	if (verbosity) {
		cout << "  data have been read from file " << fname << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

nh3_oh2_coll_data::nh3_oh2_coll_data(const string path, const energy_diagram *di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j;
	string fname;
	ifstream input;

	if (rounding(2.*di->mol.spin) == 3) 
		fname = path + "coll_nh3/coll_onh3_oh2.txt";
	else fname = path + "coll_nh3/coll_pnh3_oh2.txt";

	input.open(fname.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lev >> jmax;
	jmax++;      // one point is reserved for 0 K;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));
	
	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	for (i = 0; i < imax; i++)
	{
		input >> j >> j >> j;
		for (j = 1; j < jmax; j++) {
			input >> coeff[i][j];
		}
	}
	input.close();
	calc_coeff_deriv();
	
	if (verbosity) {
		cout << "  data have been read from file " << fname << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

nh3_h_coll_data::nh3_h_coll_data(const string path, const energy_diagram *di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j;
	string fname;
	ifstream input;

	if (rounding(2.*di->mol.spin) == 3) 
		fname = path + "coll_nh3/coll_onh3_h.txt";
	else fname = path + "coll_nh3/coll_pnh3_h.txt";

	input.open(fname.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lev >> jmax;
	jmax++;      // one point is reserved for 0 K;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));
	
	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	for (i = 0; i < imax; i++)
	{
		input >> j >> j >> j;
		for (j = 1; j < jmax; j++) {
			input >> coeff[i][j];
		}
	}
	input.close();
	calc_coeff_deriv();
	
	if (verbosity) {
		cout << "  data have been read from file " << fname << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}


//
// The class calculates collisional rates
//

nh3_collisions::nh3_collisions(const string &data_path, const energy_diagram* nh3_di, int verbosity)
{
	if (verbosity) 
		cout << "NH3 collisional rate coefficients are being initializing..." << endl;
	
	nb_lev = nh3_di->nb_lev;

	coll_data.push_back( new nh3_he_coll_data(data_path, nh3_di, verbosity) );
	coll_data.push_back( new nh3_ph2_coll_data(data_path, nh3_di, verbosity) );
	coll_data.push_back( new nh3_oh2_coll_data(data_path, nh3_di, verbosity) );
	coll_data.push_back( new nh3_h_coll_data(data_path, nh3_di, verbosity) );
	nb1 = (int) coll_data.size();

	// the data on electron collisions must be here;
	nb2 = (int) coll_data.size();

	// the data on H+ collisions must be here;
	nb3 = (int) coll_data.size();

	max_temp = new double [nb3];
	for (int i = 0; i < nb3; i++) {
		max_temp[i] = coll_data[i]->get_max_temp();
	}
}

void nh3_collisions::set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, 
	double h_conc, double el_conc, double *&concentration, int *&indices) const
{
	collisional_transitions::set_gas_param(temp_neutrals, temp_el, he_conc, ph2_conc, oh2_conc, h_conc, el_conc, concentration, 
		indices);

	concentration[0] = he_conc;
	concentration[1] = ph2_conc;
	concentration[2] = oh2_conc;
	concentration[3] = h_conc;
}

// the energy of the first level is higher, up_lev.nb > low_lev.nb;
void nh3_collisions::get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, 
	double &up_rate, double temp_neutrals, const double *concentration, const int *indices) const
{
	up_rate = down_rate = 0.;	
	if (up_lev.nb < coll_data[0]->nb_lev) {
		down_rate = coll_data[0]->get_rate(up_lev.nb, low_lev.nb, indices[0], (temp_neutrals < max_temp[0]) ? temp_neutrals : max_temp[0]) *concentration[0];
	}
	if (up_lev.nb < coll_data[1]->nb_lev) {
		down_rate += coll_data[1]->get_rate(up_lev.nb, low_lev.nb, indices[1], (temp_neutrals < max_temp[1]) ? temp_neutrals : max_temp[1]) *concentration[1] 
			+ coll_data[2]->get_rate(up_lev.nb, low_lev.nb, indices[2], (temp_neutrals < max_temp[2]) ? temp_neutrals : max_temp[2]) *concentration[2]
			+ coll_data[3]->get_rate(up_lev.nb, low_lev.nb, indices[3], (temp_neutrals < max_temp[3]) ? temp_neutrals : max_temp[3]) *concentration[3];
	}
	if (down_rate > MIN_COLLISION_RATE)
		up_rate = down_rate *exp((low_lev.energy - up_lev.energy)*CM_INVERSE_TO_KELVINS/temp_neutrals) *up_lev.g /((double) low_lev.g);
	else down_rate = 0.;
}
