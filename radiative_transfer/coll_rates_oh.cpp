
// 27.04.2018
// 01.04.2020 - check for errors in HF data;
//

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory.h>
#include <cmath>
#include <cfloat>

#include "coll_rates_oh.h"
#include "utils.h"
#include "constants.h"

#define MAX_TEXT_LINE_WIDTH 240
#define SOURCE_NAME "coll_rates_oh.cpp"
using namespace std;

oh_h2_coll_data::oh_h2_coll_data(const string path, const energy_diagram *di, bool is_h2_j0, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j;
	
	string fname;
	ifstream input;

	if (is_h2_j0) 
		fname = path + "coll_oh/coll_oh_h2j0.txt";
	else fname = path + "coll_oh/coll_oh_h2j1.txt";	

	input.open(fname.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lev >> jmax;
	jmax++; // one point is reserved for 0 K;

	// nb of rows in the array:
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

oh_he_coll_data::oh_he_coll_data(const string path, const energy_diagram *di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j;
	string fname;
	ifstream input;

	fname = path + "coll_oh/coll_oh_he.txt";
	input.open(fname.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lev >> jmax;
	nb_lev -= 2; // last two levels have higher numbers according to HITRAN 2016 database; they are excluded;
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

oh_hf_h2_coll_data::oh_hf_h2_coll_data(const std::string path, const energy_diagram *di, bool coll_partner_is_ortho, int verbosity)
{
    char text_line[MAX_TEXT_LINE_WIDTH];
    int i, j, li, lf, nb;
    string fname;
    ifstream input;

    if (coll_partner_is_ortho)
        fname = path + "coll_oh/coll_oh_hf_oh2.txt";
    else fname = path + "coll_oh/coll_oh_hf_ph2.txt";

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
    imax = nb_lev * (nb_lev - 1) / 2;

    tgrid = new double[jmax];
    coeff = alloc_2d_array<double>(imax, jmax);
    memset(*coeff, 0, jmax * imax * sizeof(double));

    coeff_deriv = alloc_2d_array<double>(imax, jmax);
    memset(*coeff_deriv, 0, imax * jmax * sizeof(double));

    tgrid[0] = 0.;
    for (j = 1; j < jmax; j++) {
        input >> tgrid[j];
    }

    for (i = 0; i < imax; i++)
    {
        input >> j >> li >> lf;
        nb = (li - 2) * (li - 1) / 2 + lf - 1; // the numeration of levels in the file starts from 1

        for (j = 1; j < jmax; j++) {
            input >> coeff[nb][j];
        }
    }
    input.close();
    calc_coeff_deriv();

    if (verbosity) {
        cout << "  data have been read from file " << fname << endl
            << "  temperature range " << (int)tgrid[1] << " - " << (int)tgrid[jmax - 1] << endl;
    }
}

oh_hf_he_coll_data::oh_hf_he_coll_data(const std::string path, const energy_diagram*, int verbosity)
{
    char text_line[MAX_TEXT_LINE_WIDTH];
    int i, j, li, lf, nb;
    
    string fname;
    ifstream input;

    fname = path + "coll_oh/coll_oh_hf_he.txt";
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
    imax = nb_lev * (nb_lev - 1) / 2;

    tgrid = new double[jmax];
    coeff = alloc_2d_array<double>(imax, jmax);
    memset(*coeff, 0, jmax * imax * sizeof(double));

    coeff_deriv = alloc_2d_array<double>(imax, jmax);
    memset(*coeff_deriv, 0, imax * jmax * sizeof(double));

    tgrid[0] = 0.;
    for (j = 1; j < jmax; j++) {
        input >> tgrid[j];
    }

    for (i = 0; i < imax; i++)
    {
        input >> li >> lf >> j >> j;
        nb = (li - 2) * (li - 1) / 2 + lf - 1; // the numeration of levels in the file starts from 1

        for (j = 1; j < jmax; j++) {
            input >> coeff[nb][j];
        }
    }
    input.close();
    calc_coeff_deriv();

    if (verbosity) {
        cout << "  data have been read from file " << fname << endl
            << "  temperature range " << (int)tgrid[1] << " - " << (int)tgrid[jmax - 1] << endl;
    }
}

//
// The class calculates collisional rates
//

// Without HF splitting
oh_collisions::oh_collisions(const string &data_path, const energy_diagram* oh_di, int verbosity)
{
	bool is_h2_j0;
	if (verbosity) 
		cout << "OH collisional rate coefficients are being initializing..." << endl;
	
	nb_lev = oh_di->nb_lev;

	coll_data.push_back( new oh_he_coll_data(data_path, oh_di, verbosity) );
	coll_data.push_back( new oh_h2_coll_data(data_path, oh_di, is_h2_j0 = true, verbosity) );
	coll_data.push_back( new oh_h2_coll_data(data_path, oh_di, is_h2_j0 = false, verbosity) );
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

void oh_collisions::set_gas_param(double temp_neutrals, double temp_el, double he_conc, double h2j0_conc, double h2j1_conc, 
	double h_conc, double el_conc, double *&concentration, int *&indices) const
{
	collisional_transitions::set_gas_param(temp_neutrals, temp_el, he_conc, h2j0_conc, h2j1_conc, h_conc, el_conc, concentration, 
		indices);

	concentration[0] = he_conc;
	concentration[1] = h2j0_conc;
	concentration[2] = h2j1_conc;
}

// The energy of the first level is higher, up_lev.nb > low_lev.nb;
void oh_collisions::get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, 
	double &up_rate, double temp_neutrals, const double *concentration, const int *indices) const
{
	up_rate = down_rate = 0.;
	// for oh-he collisions, 44 lower levels of oh are considered, there is no check for up.lev number,
	down_rate = coll_data[0]->get_rate(up_lev.nb, low_lev.nb, indices[0], (temp_neutrals < max_temp[0]) ? temp_neutrals : max_temp[0]) *concentration[0];
		
	if (up_lev.nb < coll_data[1]->nb_lev) {
		down_rate += coll_data[1]->get_rate(up_lev.nb, low_lev.nb, indices[1], (temp_neutrals < max_temp[1]) ? temp_neutrals : max_temp[1]) *concentration[1] 
			+ coll_data[2]->get_rate(up_lev.nb, low_lev.nb, indices[2], (temp_neutrals < max_temp[2]) ? temp_neutrals : max_temp[2]) *concentration[2];
	}
	
	if (down_rate > MIN_COLLISION_RATE)
		up_rate = down_rate *exp((low_lev.energy - up_lev.energy)*CM_INVERSE_TO_KELVINS/temp_neutrals) *up_lev.g /((double) low_lev.g);
	else down_rate = 0.;
}

// for data taking into account HF splitting
oh_hf_collisions::oh_hf_collisions(const string& data_path, const energy_diagram* oh_di, int verbosity)
{
    bool coll_partner_is_ortho;
    if (verbosity)
        cout << "OH collisional rate coefficients (including HF splitting) are being initializing..." << endl;

    nb_lev = oh_di->nb_lev;

    coll_data.push_back(new oh_hf_he_coll_data(data_path, oh_di, verbosity));
    coll_data.push_back(new oh_hf_h2_coll_data(data_path, oh_di, coll_partner_is_ortho = false, verbosity));
    coll_data.push_back(new oh_hf_h2_coll_data(data_path, oh_di, coll_partner_is_ortho = true, verbosity));
    nb1 = (int)coll_data.size();

    // the data on electron collisions must be here;
    nb2 = (int)coll_data.size();

    // the data on H+ collisions must be here;
    nb3 = (int)coll_data.size();

    max_temp = new double[nb3];
    for (int i = 0; i < nb3; i++) {
        max_temp[i] = coll_data[i]->get_max_temp();
    }
}

void oh_hf_collisions::set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc,
    double h_conc, double el_conc, double*& concentration, int*& indices) const
{
    collisional_transitions::set_gas_param(temp_neutrals, temp_el, he_conc, ph2_conc, oh2_conc, h_conc, el_conc, concentration,
        indices);

    concentration[0] = he_conc;
    concentration[1] = ph2_conc;
    concentration[2] = oh2_conc;
}

// the energy of the first level is higher, up_lev.nb > low_lev.nb;
void oh_hf_collisions::get_rate_neutrals(const energy_level& up_lev, const energy_level& low_lev, double& down_rate,
    double& up_rate, double temp_neutrals, const double* concentration, const int* indices) const
{
    up_rate = down_rate = 0.;
    if (up_lev.nb < coll_data[0]->nb_lev) {
        down_rate = coll_data[0]->get_rate(up_lev.nb, low_lev.nb, indices[0], (temp_neutrals < max_temp[0]) ? temp_neutrals : max_temp[0]) * concentration[0];
    }
    if (up_lev.nb < coll_data[1]->nb_lev) {
        down_rate += coll_data[1]->get_rate(up_lev.nb, low_lev.nb, indices[1], (temp_neutrals < max_temp[1]) ? temp_neutrals : max_temp[1]) * concentration[1]
            + coll_data[2]->get_rate(up_lev.nb, low_lev.nb, indices[2], (temp_neutrals < max_temp[2]) ? temp_neutrals : max_temp[2]) * concentration[2];
    }

    if (down_rate > MIN_COLLISION_RATE)
        up_rate = down_rate * exp((low_lev.energy - up_lev.energy) * CM_INVERSE_TO_KELVINS / temp_neutrals) * up_lev.g / ((double)low_lev.g);
    else down_rate = 0.;
}
