
//
// Modifications:
// 06.06.2018. Update is done.
//

#include <stdlib.h>
#include <memory.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <cfloat>

#include "constants.h"
#include "utils.h"
#include "coll_rates_ch3oh.h"

#define MAX_TEXT_LINE_WIDTH 240
#define SOURCE_NAME "coll_rates_ch3oh.cpp"
using namespace std;

//
// The classes contain collisional data.
//

ch3oh_he_coll_data::ch3oh_he_coll_data(const string &path, const energy_diagram *ch3oh_levels, int verbosity)
{
	char ch, text_line[MAX_TEXT_LINE_WIDTH];	
	int i, j, k, l, i1, i2, nb_of_levels, vt, vt_max;
	int *nb_arr;
	double a, rate, energy;
	
	string file_name;
	stringstream s;
	ifstream input;

	nb_lev = ch3oh_levels->nb_lev;
	imax = nb_lev*(nb_lev-1)/2;
	jmax = 41; // 40 temperature values in the CH3OH-He rovibrational collision data, +1 value for T = 0 K;

	tgrid = new double [jmax];
	tgrid[0] = 0.;

	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	coeff_deriv = alloc_2d_array<double>(imax, jmax-1);

	// The table of energy levels is listed in the file, the number of levels is 256 here;
	nb_of_levels = 256;
	nb_arr = new int [nb_of_levels];
	memset(nb_arr, 0, nb_of_levels*sizeof(int));
	
	vt_max = NB_VIBR_EXCIT_CH3OH_RABLI;
	for (vt = 0; vt <= vt_max; vt++)
	{
		s.str("");
		if (rounding(2.*ch3oh_levels->mol.spin) == 3) s << 'a';
		else s << 'e';
		s << vt;
		
		file_name = path + "coll_ch3oh/coll_ch3oh_";
		file_name += s.str();
		file_name += "_he.txt";
		input.open(file_name.c_str(), ios_base::in);

		if (!input) {
			cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
			exit(1);
		}
		
		// the four first lines are comments; 
		for (i = 0; i < 4; i++) {
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		}
		for (i = 0; i < nb_of_levels; i++)
		{
			input >> l >> ch;
			if (rounding(2.*ch3oh_levels->mol.spin) == 3) input >> ch;
			input >> j >> k >> energy;
			
			if (ch == '-') k = -k;
			nb_arr[i] = ch3oh_levels->get_nb(vt, j, k);
		}
		// the end of the previous text line is read;
		for (i = 0; i < 4; i++) {
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		}
		for (j = 1; j < 21; j++)
		{
			input >> tgrid[j];
			for (i = 0; i < nb_of_levels; i++) // row nb is an index of final level;
			{
				input >> l;
				for (k = 0; k < nb_of_levels; k++) // column nb is an index of initial level;
				{
					input >> rate;
					// some values are in the form a.b-dfg without E symbol; 
					if (fabs(rate) >= 1. - DBL_EPSILON) {
						input >> a; 
						rate = 0.;
					}
					else if (rate < 0.) 
						rate = 0.;
						
					i1 = nb_arr[k]; // initial
					i2 = nb_arr[i]; // final

					if (i1 != -1 && i2 != -1)
					{
						if (i1 > i2) 
						{
							l = i1*(i1-1)/2 + i2;
							coeff[l][j] += 0.5*rate;
						}
						else if (i1 < i2)
						{
							l = i2*(i2-1)/2 + i1;
							coeff[l][j] += 0.5*rate *ch3oh_levels->lev_array[i1].g /((double) ch3oh_levels->lev_array[i2].g)
								*exp((ch3oh_levels->lev_array[i2].energy - ch3oh_levels->lev_array[i1].energy) *CM_INVERSE_TO_KELVINS/tgrid[j]);
						}
					}
				}
			}
		}
		input.close();
		if (verbosity) {
			cout << left << "  data have been read from file " << file_name << endl 
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[20] << endl;
		}	
	}
	// The extrapolation of rate coefficients to higher temperatures;
	for (i = 0; i < imax; i++)
	{
		for (j = 21; j < jmax; j++)
			coeff[i][j] = coeff[i][20]; 
	}

	// Further there are collisional coefficients, that include ro-vibrational transitions; 
	nb_of_levels = 150; // new value of the level number;
	
	s.str("");
	if (rounding(2.*ch3oh_levels->mol.spin) == 3) s << 'a';
	else s << 'e';
		
	file_name = path + "coll_ch3oh/coll_ch3oh_";
	file_name += s.str();
	file_name += "_he_rovibr.txt";
	input.open(file_name.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	
	// four first lines are comments; 
	for (i = 0; i < 4; i++) {
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	}
	for (i = 0; i < nb_of_levels; i++)
	{
		input >> l >> vt >> ch;
		if (rounding(2.*ch3oh_levels->mol.spin) == 3) input >> ch;
		input >> j >> k >> energy;

		if (ch == '-') k = -k;
		nb_arr[i] = ch3oh_levels->get_nb(vt, j, k);
	}
	// the end of the previous text line is read;
	for (i = 0; i < 4; i++) {
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	}
	for (j = 1; j < jmax; j++)
	{
		input >> tgrid[j];
		for (i = 0; i < nb_of_levels; i++) // row nb is an index of final level;
		{
			input >> l;
			for (k = 0; k < nb_of_levels; k++) // column nb is an index of initial level;
			{
				input >> rate;					
				// some values are in the form a.b-dfg without E symbol; 
				if (fabs(rate) >= 1. - DBL_EPSILON) {
					input >> a; 
					rate = 0.;
				}
				else if (rate < 0.) 
					rate = 0.;
					
				i1 = nb_arr[k]; // initial
				i2 = nb_arr[i]; // final

				if (i1 != -1 && i2 != -1)
				{
					if (i1 > i2) 
					{
						l = i1*(i1-1)/2 + i2;
						if (k > i && ch3oh_levels->lev_array[i1].v == ch3oh_levels->lev_array[i2].v) 
							coeff[l][j] = 0.; // The data from previous files must be overriden here;
						
						coeff[l][j] += 0.5*rate;
					}
					else if (i1 < i2)
					{
						l = i2*(i2-1)/2 + i1;
						if (k > i && ch3oh_levels->lev_array[i1].v == ch3oh_levels->lev_array[i2].v) 
							coeff[l][j] = 0.;

						coeff[l][j] += 0.5*rate *ch3oh_levels->lev_array[i1].g /((double) ch3oh_levels->lev_array[i2].g)
							*exp((ch3oh_levels->lev_array[i2].energy - ch3oh_levels->lev_array[i1].energy) *CM_INVERSE_TO_KELVINS/tgrid[j]);
					}
				}
			}
		}
	}
	input.close();
	calc_coeff_deriv();

	if (verbosity) {
		cout << left << "  data have been read from file " << file_name << endl 
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
	delete [] nb_arr;
}


ch3oh_ph2_coll_data::ch3oh_ph2_coll_data(const string &path, const energy_diagram *ch3oh_levels, int verbosity)
{
	char ch, text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, k, l, i1, i2, nb_of_levels, vt, vt_max;
	int *nb_arr;
	double a, rate, energy;
	
	string file_name;
	stringstream s;
	ifstream input;

	nb_lev = ch3oh_levels->nb_lev;
	imax = nb_lev*(nb_lev-1)/2;
	jmax = 21; // one point is reserved for 0 K;

	tgrid = new double [jmax];
	tgrid[0] = 0.;

	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	coeff_deriv = alloc_2d_array<double>(imax, jmax-1);

	// The table of energy levels is listed in the file, the number of levels is 256 here;
	nb_of_levels = 256;
	nb_arr = new int [nb_of_levels];
	memset(nb_arr, 0, nb_of_levels*sizeof(int));
	
	vt_max = NB_VIBR_EXCIT_CH3OH_RABLI;
	for (vt = 0; vt <= vt_max; vt++)
	{
		s.str("");
		if (rounding(2.*ch3oh_levels->mol.spin) == 3) s << 'a';
		else s << 'e';
		s << vt;
		
		file_name = path + "coll_ch3oh/coll_ch3oh_";
		file_name += s.str();
		file_name += "_ph2.txt";
		input.open(file_name.c_str(), ios_base::in);

		if (!input) {
			cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
			exit(1);
		}
		
		// four first lines are comments;
		for (i = 0; i < 4; i++) {
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		}
		for (i = 0; i < nb_of_levels; i++)
		{
			input >> l >> ch;
			if (rounding(2.*ch3oh_levels->mol.spin) == 3) input >> ch;
			input >> j >> k >> energy;

			if (ch == '-') k = -k;
			nb_arr[i] = ch3oh_levels->get_nb(vt, j, k);
		}
		// the end of the previous text line is read;
		for (i = 0; i < 4; i++) {
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		}
		for (j = 1; j < jmax; j++)
		{
			input >> tgrid[j];
			for (i = 0; i < nb_of_levels; i++) // row nb is an index of final level;
			{
				input >> l;
				for (k = 0; k < nb_of_levels; k++) // column nb is an index of initial level;
				{
					input >> rate; 
					// some values are in the form a.b-dfg without E symbol; 
					if (fabs(rate) >= 1.-DBL_EPSILON) {
						input >> a; 
						rate = 0.;
					}
					else if (rate < 0.) 
						rate = 0.;
						
					i1 = nb_arr[k]; // initial
					i2 = nb_arr[i]; // final

					if (i1 != -1 && i2 != -1)
					{
						if (i1 > i2) 
						{
							l = i1*(i1-1)/2 + i2;
							coeff[l][j] += 0.5*rate;
						}
						else if (i1 < i2)
						{
							l = i2*(i2-1)/2 + i1;
							coeff[l][j] += 0.5*rate *ch3oh_levels->lev_array[i1].g /((double) ch3oh_levels->lev_array[i2].g)
								*exp((ch3oh_levels->lev_array[i2].energy - ch3oh_levels->lev_array[i1].energy) *CM_INVERSE_TO_KELVINS/tgrid[j]);
						}
					}
				}
			}
		}
		input.close();
		
		if (verbosity) {
			cout << left << "  data have been read from file " << file_name << endl 
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
		}
	}
	calc_coeff_deriv();
	delete [] nb_arr;
}


ch3oh_oh2_coll_data::ch3oh_oh2_coll_data(const string &path, const energy_diagram *ch3oh_levels, int verbosity)
{
	char ch, text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, k, l, i1, i2, nb_of_levels;
	int *nb_arr;
	double a, rate, energy;
	
	string file_name;
	stringstream s;
	ifstream input;

	nb_lev = ch3oh_levels->nb_lev;
	imax = nb_lev*(nb_lev-1)/2;
	jmax = 21; // one point is reserved for 0 K;

	tgrid = new double [jmax];
	tgrid[0] = 0.;

	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	coeff_deriv = alloc_2d_array<double>(imax, jmax-1);

	// The table of energy levels is listed in the file, the number of levels is 100 for these data;
	nb_of_levels = 100;
	nb_arr = new int [nb_of_levels];
	memset(nb_arr, 0, nb_of_levels*sizeof(int));
	
	s.str("");
	if (rounding(2.*ch3oh_levels->mol.spin) == 3) s << "a0";
	else s << "e0";
	
	file_name = path + "coll_ch3oh/coll_ch3oh_";
	file_name += s.str();
	file_name += "_oh2.txt";
	input.open(file_name.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
		
	// four first lines are comments;
	for (i = 0; i < 4; i++) {
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	}
	for (i = 0; i < nb_of_levels; i++)
	{
		input >> l >> ch;
		if (rounding(2.*ch3oh_levels->mol.spin) == 3) input >> ch;
		input >> j >> k >> energy;

		if (ch == '-') k = -k;
		nb_arr[i] = ch3oh_levels->get_nb(0, j, k);
	}
	// the end of the previous text line is read;
	for (i = 0; i < 4; i++) {
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	}
	for (j = 1; j < jmax; j++)
	{
		input >> tgrid[j];
		for (i = 0; i < nb_of_levels; i++) // row nb is an index of final level;
		{
			input >> l;
			for (k = 0; k < nb_of_levels; k++) // column nb is an index of initial level;
			{
				input >> rate;
				// some values are in the form a.b-dfg without E symbol; 
				if (fabs(rate) >= 1. - DBL_EPSILON) {
					input >> a; 
					rate = 0.;
				}
				else if (rate < 0.) 
					rate = 0.;
						
				i1 = nb_arr[k]; // initial
				i2 = nb_arr[i]; // final

				if (i1 != -1 && i2 != -1)
				{
					if (i1 > i2) 
					{
						l = i1*(i1-1)/2 + i2;
						coeff[l][j] += 0.5*rate;
					}
					else if (i1 < i2)
					{
						l = i2*(i2-1)/2 + i1;
						coeff[l][j] += 0.5*rate *ch3oh_levels->lev_array[i1].g /((double) ch3oh_levels->lev_array[i2].g)
							*exp((ch3oh_levels->lev_array[i2].energy - ch3oh_levels->lev_array[i1].energy) *CM_INVERSE_TO_KELVINS/tgrid[j]);
					}
				}
			}
		}
	}
	input.close();
	calc_coeff_deriv();

	if (verbosity) {
		cout << left << "  data have been read from file " << file_name << endl 
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
	delete [] nb_arr;
}

//
// The class that calculates the collisional rates of CH3OH
//

ch3oh_collisions::ch3oh_collisions(const string &path, const energy_diagram *ch3oh_levels, int verbosity) 
{
	if (verbosity) 
		cout << "CH3OH collisional rate coefficients are being initializing..." << endl;
	
	nb_lev = ch3oh_levels->nb_lev;

	// data for rovibrational transitions includes temperatures >= 200 K (for H2) and >= 120 K (for He);
	coll_data.push_back( new ch3oh_he_coll_data(path, ch3oh_levels, verbosity) );
	coll_data.push_back( new ch3oh_ph2_coll_data(path, ch3oh_levels, verbosity) );
	coll_data.push_back( new ch3oh_oh2_coll_data(path, ch3oh_levels, verbosity) );
	nb1 = (int) coll_data.size();

	// data on electron collisions must be here;
	nb2 = (int) coll_data.size();

	// data on H+ collisions must be here;
	nb3 = (int) coll_data.size();

	max_temp = new double [nb3];
	for (int i = 0; i < nb3; i++) {
		max_temp[i] = coll_data[i]->get_max_temp();
	}
}

void ch3oh_collisions::set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, double h_conc,
		double el_conc, double *&concentration, int *&indices) const
{
	collisional_transitions::set_gas_param(temp_neutrals, temp_el, he_conc, ph2_conc, oh2_conc, h_conc, el_conc, concentration, 
		indices);

	concentration[0] = he_conc;
	concentration[1] = ph2_conc;
	concentration[2] = oh2_conc;
}

// The energy of the first level is higher, up_lev.nb > low_lev.nb
void ch3oh_collisions::get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, double &up_rate,
		double temp_neutrals, const double *concentration, const int *indices) const
{
	// The level numbers are in the range of the data arrays by default, there is no upper restriction on the level nb;
	// in the CH3OH-H2 collision data only the rates for torsionally elastic transitions are present;
	// for ortho-H2 the rates of transitions including levels v_t = 0 and J <= 9 are available;  
	if (up_lev.v == low_lev.v) {
		if (up_lev.v == 0 && up_lev.j <= 9 && low_lev.j <= 9) // in this case the data for ortho-H2 are available;
		{
			down_rate = coll_data[1]->get_rate(up_lev.nb, low_lev.nb, indices[1], (temp_neutrals < max_temp[1]) ? temp_neutrals : max_temp[1]) *concentration[1] 
				+ coll_data[2]->get_rate(up_lev.nb, low_lev.nb, indices[2], (temp_neutrals < max_temp[2]) ? temp_neutrals : max_temp[2]) *concentration[2];
		}
		else
		{ // in this case the data for ortho-H2 are assumed to be identical to para-H2;
			down_rate = coll_data[1]->get_rate(up_lev.nb, low_lev.nb, indices[1], (temp_neutrals < max_temp[1]) ? temp_neutrals : max_temp[1]) 
				*(concentration[1] + concentration[2]);
		}
		// CH3OH-He
		down_rate += coll_data[0]->get_rate(up_lev.nb, low_lev.nb, indices[0], (temp_neutrals < max_temp[0]) ? temp_neutrals : max_temp[0]) *concentration[0];
	}
	else {
	// For torsionally inelastic transition induced by para-H2, Rabli, Flower (2011) suggested using the rate coefficients for He, 
	// and three times these values for collisions with ortho-H2. 
	// initial data sets are available for <= 400 K;
		down_rate = coll_data[0]->get_rate(up_lev.nb, low_lev.nb, indices[0], (temp_neutrals < max_temp[0]) ? temp_neutrals : max_temp[0]) 
				*(concentration[0] + concentration[1] + 3.*concentration[2]);
	}

	if (down_rate > MIN_COLLISION_RATE)
		up_rate = down_rate *exp((low_lev.energy - up_lev.energy)*CM_INVERSE_TO_KELVINS/temp_neutrals) *up_lev.g /((double) low_lev.g);
	else up_rate = down_rate = 0.;
}
