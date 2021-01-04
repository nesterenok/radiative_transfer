// The version #2 of the spectroscopy source file.
/* Modifications:
	13.05.2016. The error was found in the constructor of h2_diagram and h2o_diagram: the variable nb in the constructor 
		was introduced, i++ was moved outside the if condition;
	16.05.2016. Comment lines were added to the file levels_h2.txt;
	01.06.2016. The routines save_populations and read_populations were added (were copied from the iteration_control.cpp file); 
	17.06.2016. The check was added in the constructor of the class for H2 Einstein coefficients; 
	20.06.2016. The classes for atoms and ions were added; 
	21.06.2016. The parameter "name" was added to the class energy_level;
	07.07.2016. The h2_diagram* was changed by the energy_diagram* in the constructor variable lists;
	19.09.2016. The class energy_diagram now contains an object molecule, not the pointer to this object;
	28.10.2016. The class for CO molecule was written;
	03.03.2017. The error was found in methanol energy level initialisation, minimal energy was set erroneous.
		Check for errors.
	08.09.2017. Two levels were added to the H2 level list, that are present in the Dabrowski (1984); 
		v=11 j=13 e=36310.61 cm-1; v=11 j=14 e=36593.4 cm-1. Spectroscopic data on CO were updated, HITRAN2016 is used.
		Check for errors.
	25.04.2018. Angular momentum and its projections on molecular axes are float pointing values. OH molecule data were added.
	10.05.2018. NH3 molecule data were added.
    06.04.2020. The energy_level class, operator == was changed.
	09.06.2020. The statistical weight calculation is done unique in all classes (the j and spin are taken into account).
*/

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <stdlib.h>
#include <memory.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <cfloat>

#include "constants.h"
#include "utils.h"
#include "spectroscopy.h"

#define MAX_TEXT_LINE_WIDTH 240  // maximal size of the comment lines in the files,
#define SOURCE_NAME "spectroscopy.cpp"
using namespace std;

molecule::molecule(const string &n, int is, double m, double sp) 
: name(n), isotop(is), mass(m), spin(sp)
{;}

//
// The classes that contain data on energy levels.
//

energy_level::energy_level()
: nb(0), v(0), syminv(0), g(0), j(0.), k1(0.), k2(0.), spin(0.), energy(0.), hf(0.), name("")
{;}

energy_diagram::energy_diagram(molecule m, int verb) : mol(m), nb_lev(0), verbosity(verb), hyperfine_splitting(false)
{;}

energy_diagram::~energy_diagram() {
	lev_array.clear();
}

void energy_diagram::report(const std::string & fname)
{
    if (verbosity > 0) {
        cout << "Specimen levels have been initialized: " << mol.name << endl
            << "    data have been read from file: " << fname << endl
            << "    number of levels: " << nb_lev << endl;
    }
}

h2_diagram::h2_diagram(const string &data_path, molecule m, int &n_l, int verb) : energy_diagram(m, verb)
{	
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, i_max, v, j, nb;
	double energy;

	string fname;
	ifstream input;
	energy_level level;

	fname = data_path + "spectroscopy/levels_h2.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}
	
	// Comment lines:
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	nb = i = 0;
	input >> i_max;

	while (nb < n_l && i < i_max)
	{
		input >> v >> j >> energy;
		if (j%2 == rounding(mol.spin) || rounding(mol.spin) == -1)
		{
			level.v = v;
			level.j = j;
			level.energy = energy;
				
			level.spin = j%2; // mol.spin value can be undefined
			level.g = (2*j + 1) *(2*(j%2) + 1);
			level.nb = nb;
			
			lev_array.push_back(level);
			nb++;
		}
		i++;
	}
	input.close();
	nb_lev = n_l = (int) lev_array.size();
    report(fname);
}

int h2_diagram::get_nb(int v, double j) const
{
	for (int i = 0; i < nb_lev; i++) {
		if (lev_array[i].v == v && rounding(lev_array[i].j) == rounding(j))
			return i;
	}
	return -1;
}

h2o_diagram::h2o_diagram(const string &data_path, molecule m, int &n_l, int nb_vibr, int verb) : energy_diagram(m, verb)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, i_max, v, v1, v2, v3, j, ka, kc, nb;	
	double energy;

	string fname;
	ifstream input;
	energy_level level;
	
	if (nb_vibr > NB_VIBR_EXCIT_H2O) 
		nb_vibr = NB_VIBR_EXCIT_H2O;
	
	if (mol.isotop == 1) 
		fname = data_path + "spectroscopy/levels_h2o16.txt";
	else fname = data_path + "spectroscopy/levels_h2o18.txt";

	input.open(fname.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> i_max;

	nb = i = 0;
	while (nb < n_l && i < i_max)
	{
		input >> v1 >> v2 >> v3;
		input >> j >> ka >> kc >> energy;  //j, ka, kc are integers
			
		v = get_vibr_nb(v1, v2, v3);
		if (abs(ka + kc + v3)%2 == rounding(mol.spin) && v <= nb_vibr)
		{ 	
			level.v = v;			
			level.j = j;
			level.k1 = ka;
			level.k2 = kc;
			level.energy = energy;
			level.spin = mol.spin;

			// either ortho- or para-H2O molecules is considered, 
			level.g = rounding(2. * mol.spin + 1.) * (2 * j + 1);
			level.nb = nb;

			lev_array.push_back(level);
			nb++;
		}
		i++;
	}
	input.close();
	nb_lev = n_l = (int) lev_array.size();
    report(fname);
}

int h2o_diagram::get_nb(int v, double j, double tau) const
{
	for (int i = 0; i < nb_lev; i++) {
		if (lev_array[i].v == v && rounding(lev_array[i].j) == rounding(j) && rounding(lev_array[i].k1-lev_array[i].k2) == rounding(tau))
			return i;
	}
	return -1;
}

int h2o_diagram::get_vibr_nb(int v1, int v2, int v3) const
{
	if (v1 == 0 && v2 == 0 && v3 == 0) return 0;
	else if (v1 == 0 && v2 == 1 && v3 == 0) return 1;
	else if (v1 == 0 && v2 == 2 && v3 == 0) return 2;
	else if (v1 == 1 && v2 == 0 && v3 == 0)	return 3;
	else if (v1 == 0 && v2 == 0 && v3 == 1)	return 4;

	return 5;
}

ch3oh_diagram::ch3oh_diagram(const string &data_path, molecule m, int &n_l, int nb_vibr, int ang_mom_max, int verb) 
	: energy_diagram(m, verb)
{
	char ch1, ch2, text_line[MAX_TEXT_LINE_WIDTH];
	int i, l, j, k1, vt, vt_max;
	double energy, energy_min;

	string fname;
	ifstream input;
	energy_level level;
	
	fname = data_path + "spectroscopy/levels_ch3oh.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}
	
	if (nb_vibr > NB_VIBR_EXCIT_CH3OH_MEKHTIEV) 
		nb_vibr = NB_VIBR_EXCIT_CH3OH_MEKHTIEV;
		
	// two first lines are comments; 
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	// for the data given in the file - torsion-vibrational quantum number vt: 
	vt_max = 8;
		
	i = 0;
	energy_min = 0.;
	for (j = 0; j <= NB_ANG_MOM_CH3OH; j++)
	{
		for (l = 0; l < 4; l++)
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
			
		for (l = 0; l < 2*(2*j+1); l++)
		{
			input >> ch1;
			if (l < 2*j+1) input >> ch2;
			else ch2 = ' ';
				
			input >> j >> k1;
			for (vt = 0; vt <= vt_max; vt++)
			{
				input >> energy;
				if (j == 0 && k1 == 0 && vt == 0 && l == 0)
					energy_min = energy;
					
				if (vt <= nb_vibr && j <= ang_mom_max)
				{
					if ((rounding(2.*mol.spin) == 3 && ch1 == 'A') || (rounding(2.*mol.spin) == 1 && ch1 == 'E')) 
					{
						level.v = vt;
						level.j = j;
						
			// the A-symmetry species have close pairing of levels, they differ in the sign of k1 here: 
						if (ch2 == '-') 
							level.k1 = -k1;
						else level.k1 = k1;
							
						level.spin = mol.spin;
						level.energy = energy - energy_min;

			// either A- or E- CH3OH spin isomers is considered, 
						level.g = rounding(2. * mol.spin + 1.) * (2 * j + 1);
						lev_array.push_back(level);
						i++;
					}
				}
			}
		}
		input.getline(text_line, MAX_TEXT_LINE_WIDTH); // the end of the line is read;
	}
	input.close();
	
	sort(lev_array.begin(), lev_array.end());
	
	// the nb of levels can be higher or lower than the value given;
	l = (int) lev_array.size() - n_l;
	for (i = 0; i < l; i++) {
		lev_array.pop_back();
	}

	// nb of levels is corrected:
	nb_lev = n_l = (int) lev_array.size();
	for (i = 0; i < nb_lev; i++) {
		lev_array[i].nb = i;
	}
    report(fname);
}

int ch3oh_diagram::get_nb(int v, double j, double k) const
{
	for (int i = 0; i < nb_lev; i++) {
		if (lev_array[i].v == v && rounding(lev_array[i].j) == rounding(j) && rounding(lev_array[i].k1) == rounding(k))
			return i;
	}
	return -1;
}

ion_diagram::ion_diagram(const std::string &path, molecule m, int &n_l, int verb) : energy_diagram(m, verb)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i;

	string str, fname;
	ifstream input;
	energy_level level;

	fname = path + "spectroscopy/levels_";
	fname += mol.name;
	fname += ".txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	// three first lines are comments; 
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> i;
	nb_lev = n_l = (i < n_l) ? i : n_l;
	
	for (i = 0; i < nb_lev; i++)
	{
		input >> level.spin >> level.energy >> level.name >> str;
		level.name += str;
		level.nb = i;
		level.g = rounding(2.*level.spin + 1.);
		lev_array.push_back(level);
	}
	input.close();
    report(fname);
}

int ion_diagram::get_nb(const string name, int stat_w) const
{
	for (int i = 0; i < nb_lev; i++) {
		if (lev_array[i].name == name && stat_w == lev_array[i].g)
			return i;
	}
	return -1;
}

co_diagram::co_diagram(const string &path, molecule m, int &n_l, int nb_vibr, int verb) : energy_diagram(m, verb)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, i_max, v, j, nb;	
	double energy;

	string fname;
	ifstream input;
	energy_level level;
	
	if (nb_vibr > NB_VIBR_EXCIT_CO) 
		nb_vibr = NB_VIBR_EXCIT_CO;
	
	fname = path + "spectroscopy/levels_co.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> i_max;

	nb = i = 0;
	while (nb < n_l && i < i_max)
	{
		input >> v >> j >> energy;
		if (v <= nb_vibr)
		{ 	
			level.v = v;			
			level.j = j;
			level.energy = energy;
			level.spin = 0;

			level.g = 2 * j + 1;
			level.nb = nb;

			lev_array.push_back(level);
			nb++;
		}
		i++;
	}
	input.close();
	nb_lev = n_l = (int) lev_array.size();
    report(fname);
}

int co_diagram::get_nb(int v, double j) const
{
	for (int i = 0; i < nb_lev; i++) {
		if (lev_array[i].v == v && rounding(lev_array[i].j) == rounding(j))
			return i;
	}
	return -1;
}

// OH
oh_diagram::oh_diagram(const std::string &path, molecule m, int &n_l, int verb) : energy_diagram(m, verb)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i_max, v, nb, parity; 
	double j, omega, energy;

	string fname;
	ifstream input;
	energy_level level;
	
	fname = path + "spectroscopy/levels_oh.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
    input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> i_max;

	nb = 0;
	while (nb < n_l && nb < i_max)
	{
		input >> v >> j >> omega >> parity >> energy;  // j, omega are double

		level.v = v;			
		level.j = j; 
		level.k1 = omega;
		level.syminv = parity;
		level.energy = energy;
        level.spin = mol.spin;

		level.g = 2 * rounding(2.*j + 1.); // factor 2 is due to hyperfine splitting;
		level.nb = nb;

		lev_array.push_back(level);
		nb++;
	}
	input.close();
	nb_lev = n_l = (int) lev_array.size();
    report(fname);
}

// note: the angular momentum is rational number,
int oh_diagram::get_nb(int parity, int v, double j, double omega) const
{
	for (int i = 0; i < nb_lev; i++) {
		if (lev_array[i].v == v && rounding(2.*lev_array[i].j) == rounding(2.*j) && rounding(2.*lev_array[i].k1) == rounding(2.*omega)
			&& lev_array[i].syminv == parity)
			return i;
	}
	return -1;
}

oh_hf_diagram::oh_hf_diagram(const std::string& path, molecule m, int& n_l, int verb) : energy_diagram(m, verb)
{
    char text_line[MAX_TEXT_LINE_WIDTH];
    int i_max, v, nb, parity, hf;
    double j, omega, energy; // j is double here, hf is the total angular momentum including spin F = J+S (is integer)

    string fname;
    ifstream input;
    energy_level level;

    fname = path + "spectroscopy/levels_oh_hf.txt";
    input.open(fname.c_str(), ios_base::in);

    if (!input) {
        cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
        exit(1);
    }

    input.getline(text_line, MAX_TEXT_LINE_WIDTH);
    input.getline(text_line, MAX_TEXT_LINE_WIDTH);
    input.getline(text_line, MAX_TEXT_LINE_WIDTH);
    input >> i_max;

    nb = 0;
    while (nb < n_l && nb < i_max)
    {
        input >> v >> j >> omega >> parity >> hf >> energy;  // j, omega are double, hf is integer,

        level.v = v;
        level.j = j;
        level.k1 = omega;
        level.syminv = parity;
        level.hf = hf;
        level.energy = energy;
        level.spin = mol.spin;

        level.g = 2 * hf + 1; // there is no any factor;
        level.nb = nb;

        lev_array.push_back(level);
        nb++;
    }
    input.close();

    hyperfine_splitting = true;
    nb_lev = n_l = (int)lev_array.size();
    report(fname);
}

// note: the angular momentum is rational number,
int oh_hf_diagram::get_nb(int parity, int v, double j, double omega, double hf) const
{
    for (int i = 0; i < nb_lev; i++) {
        if (lev_array[i].v == v && rounding(2. * lev_array[i].j) == rounding(2. * j) && rounding(2. * lev_array[i].k1) == rounding(2. * omega)
            && rounding(2.*lev_array[i].hf) == rounding(2.*hf) && lev_array[i].syminv == parity)
            return i;
    }
    return -1;
}

// NH3
nh3_diagram::nh3_diagram(const std::string &path, molecule m, int &n_l, int verb) : energy_diagram(m, verb)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, i_max, v, nb, j, k, syminv; 
	double energy; 

	string fname;
	ifstream input;
	energy_level level;
	
	fname = path + "spectroscopy/levels_nh3.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> i_max;

	nb = i = 0;
	while (nb < n_l && i < i_max)
	{
		input >> v >> j >> k >> syminv >> energy;  // v, j, k, syminv are integers
		if ((k%3 == 0 && rounding(2.*mol.spin) == 3)  // ortho- or para-
			|| (k%3 != 0 && rounding(2.*mol.spin) == 1))
		{
			level.v = v;			
			level.j = j; 
			level.k1 = k;
			level.syminv = syminv;
			level.energy = energy;
			level.spin = mol.spin;

			// HITRAN for N^{14}H3: state independent weight 3, state dependent 4 (ortho-), 2 (para-)
			level.g = 3 * rounding(2.*mol.spin + 1.) * (2 * j + 1); // ortho = 12*(2j+1), para = 6*(2j+1)
			level.nb = nb;

			lev_array.push_back(level);
			nb++;	
		}
		i++;
	}
	input.close();
	nb_lev = n_l = (int) lev_array.size();
    report(fname);
}

int nh3_diagram::get_nb(int syminv, int v, double j, double k) const
{
	for (int i = 0; i < nb_lev; i++) {
		if (lev_array[i].v == v && rounding(lev_array[i].j) == rounding(j) && rounding(lev_array[i].k1) == rounding(k)
			&& lev_array[i].syminv == syminv)
			return i;
	}
	return -1;
}

h2co_diagram::h2co_diagram(const std::string& path, molecule m, int& n_l, int verb) : energy_diagram(m, verb)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, i_max, v, j, ka, kc, nb;
	double energy;

	string fname;
	ifstream input;
	energy_level level;

	fname = path + "spectroscopy/levels_h2co.txt";
	input.open(fname.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> i_max;

	nb = i = 0;
	while (nb < n_l && i < i_max) {
		input >> v >> j >> ka >> kc >> energy;  // v, j, ka, kc are integers

		if ((ka % 2 == 0 && rounding(mol.spin) == 0) // para- or ortho-
			|| (ka % 2 != 0 && rounding(mol.spin) == 1))
		{
			level.v = v;
			level.j = j;
			level.k1 = ka;
			level.k2 = kc;
			level.energy = energy;
			level.spin = mol.spin;

			level.g = rounding(2. * mol.spin + 1.) * (2 * j + 1); 
			level.nb = nb;

			lev_array.push_back(level);
			nb++;	
		}
		i++;
	}
	input.close();
	nb_lev = n_l = (int)lev_array.size();
	report(fname);
}

// the same as for H2O: k1 - k2 is considered as a third parameter,
int h2co_diagram::get_nb(int v, double j, double k) const
{
	for (int i = 0; i < nb_lev; i++) {
		if (lev_array[i].v == v && rounding(lev_array[i].j) == rounding(j) && rounding(lev_array[i].k1 - lev_array[i].k2) == rounding(k))
			return i;
	}
	return -1;
}

//
//	The classes that contain Einstein coefficients
//

einstein_coeff::einstein_coeff(const energy_diagram *di) : diagram(di)
{
	nb_lev = diagram->nb_lev;	
	arr = alloc_2d_array<double>(nb_lev, nb_lev);
	memset(*arr, 0, nb_lev*nb_lev*sizeof(double));
}

einstein_coeff::~einstein_coeff() {
	free_2d_array<double>(arr);
}

h2_einstein_coeff::h2_einstein_coeff(const string &path, const energy_diagram* h2_di, int verbosity) : einstein_coeff(h2_di)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int	i, i_max, v, j, up, low;
	int *check;
	double coeff, energy;

	string file_name;
	ifstream input;
	
	if (verbosity) 
		cout << "H2 molecule radiative coefficients are being initializing..." << endl;

	file_name = path + "spectroscopy/radiative_h2.txt";	
	input.open(file_name.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << file_name << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> i_max;

	check = new int [nb_lev];
	memset(check, 0, nb_lev*sizeof(int));

	for (i = 0; i < i_max; i++)
	{
		input >> v >> j;
		up = h2_di->get_nb(v, j);
			
		input >> v >> j;
		low = h2_di->get_nb(v, j);
			
		input >> coeff >> energy;
		if (low != -1 && up != -1)
		{
			arr[up][low] = coeff;
			arr[low][up] = h2_di->lev_array[up].g *coeff / ((double) h2_di->lev_array[low].g);
			check[up] = 1;
		}
	}
	input.close();
	
	if (verbosity) 
	{
		cout << "  data have been read from file " << file_name << endl;
		for (i = 0; i < nb_lev; i++) 
		{
			if (check[i] == 0) {
				cout << left << "H2 level has not radiative decay coefficient, (v,j): " << setw(5) << h2_di->lev_array[i].v 
					<< setw(5) << h2_di->lev_array[i].j << endl;
			}
		}
	}
	delete [] check;
}

h2o_einstein_coeff::h2o_einstein_coeff(const string &path, const h2o_diagram* h2o_di, int verbosity) : einstein_coeff(h2o_di)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int	i, i_max, v, v1, v2, v3, j, ka, kc, up, low;
	double coeff, energy;
	
	string file_name;
	ifstream input;
	
	if (verbosity) 
		cout << "H2O molecule radiative coefficients are being initializing..." << endl;

	if (h2o_di->mol.isotop == 1) 
		file_name = path + "spectroscopy/radiative_h2o16.txt";
	else file_name = path + "spectroscopy/radiative_h2o18.txt";

	input.open(file_name.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << file_name << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> i_max;
	
	for (i = 0; i < i_max; i++)
	{
		input >> v1 >> v2 >> v3 >> j >> ka >> kc;

		v = h2o_di->get_vibr_nb(v1, v2, v3);
		up = h2o_di->get_nb(v, j, ka-kc);
			
		input >> v1 >> v2 >> v3 >> j >> ka >> kc;
		v = h2o_di->get_vibr_nb(v1, v2, v3);
		low = h2o_di->get_nb(v, j, ka-kc);
			
		input >> coeff >> energy;
		if (low != -1 && up != -1)
		{
			arr[up][low] = coeff;
			arr[low][up] = h2o_di->lev_array[up].g *coeff / ((double) h2o_di->lev_array[low].g);
		}
	}
	input.close();
	if (verbosity) 
		cout << "  data are read from file " << file_name << endl;
}

ch3oh_einstein_coeff::ch3oh_einstein_coeff(const string &path, const energy_diagram *ch3oh_di, int verbosity) : einstein_coeff(ch3oh_di)
{
	char ch, text_line[MAX_TEXT_LINE_WIDTH];
	int i, l, j, v, k, i_max, up, low;
	double a, energy, line_strength;

	string file_name;
	ifstream input;

	if (verbosity) 
		cout << "CH3OH molecule radiative coefficients are being initializing..." << endl;
	
	if (rounding(2.*ch3oh_di->mol.spin) == 3) 
		file_name = path + "spectroscopy/radiative_ch3oh_a.txt";
	else file_name = path + "spectroscopy/radiative_ch3oh_e.txt";

	input.open(file_name.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << file_name << endl;
		exit(1);
	}
	
	// the four first lines are comments; 
	for (l = 0; l < 4; l++) {
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	}
	input >> i_max;
	for (i = 0; i < i_max; i++)
	{	
		input.get(ch); // the symbol of the line end is read;
		input.get(ch);

		input >> v >> j >> k;
		input.get(ch); // for A-species ch = +-, for E-species ch = ' ';
			
		if (ch == '-') k = -k; 
		up = ch3oh_di->get_nb(v, j, k);

		input >> v >> j >> k;
		input.get(ch);

		if (ch == '-') k = -k;
		low = ch3oh_di->get_nb(v, j, k);
		
		input >> energy >> a >> line_strength;
		for (l = 0; l < 6; l++) {
			input >> a;
		}
		if (low != -1 && up != -1) {
			energy = ch3oh_di->lev_array[up].energy - ch3oh_di->lev_array[low].energy;

			// Thorne A., Spectrophysics, 1988, p. 290-291: A (s^{-1}) = 2.819e+73 *S(C^2 m^2)/(g_up *(l(nm))^3) (?)
			// S(C^2 m^2) = (10 light_speed)^{-2} *S(CGS), l(nm)^{-3} = 10^{-21} *l(cm)^{-3}
			arr[up][low] = 64.*line_strength *DEBYE*DEBYE *M_PI *pow(M_PI*energy, 3.) 
				/(3.*PLANCK_CONSTANT *(2.*ch3oh_di->lev_array[up].j+1.));
				
			arr[low][up] = ch3oh_di->lev_array[up].g *arr[up][low]/ ((double) ch3oh_di->lev_array[low].g);
		}
	}
	input.close();
	if (verbosity) 
		cout << "  data have been read from file " << file_name << endl;
}

ion_einstein_coeff::ion_einstein_coeff(const string &path, const energy_diagram *di, int verbosity) : einstein_coeff(di)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, nb, low, up;
	double a, spin_up, spin_l;

	string name_l, name_up, str, file_name;
	stringstream ss; 
	ifstream input;

	if (verbosity) 
		cout << di->mol.name << " ion radiative coefficients are being initializing..." << endl;
	
	file_name = path + "spectroscopy/radiative_";
	file_name += di->mol.name;
	file_name += ".txt";
	input.open(file_name.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << file_name << endl;
		exit(1);
	}
	
	// three first lines are comments; 
	for (i = 0; i < 3; i++) {
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	}
	input >> nb;
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	for (i = 0; i < nb; i++)
	{
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		ss.str(text_line);

		ss >> a;
		ss.getline(text_line, MAX_TEXT_LINE_WIDTH, '|');
		ss.getline(text_line, MAX_TEXT_LINE_WIDTH, '|');
		ss.getline(text_line, MAX_TEXT_LINE_WIDTH, '|');
		
		ss >> name_l;
		ss.getline(text_line, MAX_TEXT_LINE_WIDTH, '|');
		ss >> str;
		name_l += str;
		ss.getline(text_line, MAX_TEXT_LINE_WIDTH, '|');
		
		ss >> spin_l;
		ss.getline(text_line, MAX_TEXT_LINE_WIDTH, '|');

		ss >> name_up;
		ss.getline(text_line, MAX_TEXT_LINE_WIDTH, '|');
		ss >> str;
		name_up += str;
		ss.getline(text_line, MAX_TEXT_LINE_WIDTH, '|');
		
		ss >> spin_up;

		low = di->get_nb(name_l, rounding(2.*spin_l + 1.));
		up = di->get_nb(name_up, rounding(2.*spin_up + 1.));

		if (low != -1 && up != -1) 
		{
			// The data may contain coefficients for different transition types (forbidden magnetic and electric) for the same pair of levels, 
			// these coefficients are summed;
			arr[up][low] += a;		
			arr[low][up] += di->lev_array[up].g *arr[up][low] / ((double) di->lev_array[low].g);
		}

		if (ss.eof() || ss.fail() || ss.bad())
			cout << "Error ocurred in " << SOURCE_NAME << " while reading the file " << file_name << endl;
	}
	if (verbosity) 
		cout << "  data have been read from file " << file_name << endl;
}

co_einstein_coeff::co_einstein_coeff(const string &path, const energy_diagram *di, int verbosity) : einstein_coeff(di)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int	i, i_max, v, j, up, low;
	double coeff, energy;
	
	string file_name;
	ifstream input;
	
	if (verbosity) 
		cout << "CO molecule radiative coefficients are being initializing..." << endl;
	
	file_name = path + "spectroscopy/radiative_co.txt";
	input.open(file_name.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << file_name << endl;
		exit(1);
	}	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> i_max;
	
	for (i = 0; i < i_max; i++)
	{
		input >> v >> j;
		up = di->get_nb(v, j);
			
		input >> v >> j;
		low = di->get_nb(v, j);
			
		input >> coeff >> energy;
		if (low != -1 && up != -1)
		{
			arr[up][low] = coeff;
			arr[low][up] = di->lev_array[up].g *coeff / ((double) di->lev_array[low].g);
		}
	}
	input.close();
	if (verbosity) 
		cout << "  data are read from file " << file_name << endl;
}

oh_einstein_coeff::oh_einstein_coeff(const string &path, const energy_diagram *di, int verbosity) : einstein_coeff(di)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int	i, i_max, parity, v, up, low;
	double coeff, energy, j, omega;
	
	string file_name;
	ifstream input;
	
	if (verbosity) 
		cout << "OH molecule radiative coefficients are being initializing..." << endl;
	
	file_name = path + "spectroscopy/radiative_oh.txt";
	input.open(file_name.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << file_name << endl;
		exit(1);
	}	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> i_max;
	
	for (i = 0; i < i_max; i++)
	{
		input >> v >> j >> omega >> parity;
		up = di->get_nb(parity, v, j, omega);
			
		input >> v >> j >> omega >> parity;
		low = di->get_nb(parity, v, j, omega);
			
		input >> coeff >> energy;
		if (low != -1 && up != -1)
		{
			arr[up][low] = coeff;
			arr[low][up] = di->lev_array[up].g *coeff / ((double) di->lev_array[low].g);
		}
	}
	input.close();
	if (verbosity) 
		cout << "  data are read from file " << file_name << endl;
}

oh_hf_einstein_coeff::oh_hf_einstein_coeff(const std::string& path, const energy_diagram *di, int verbosity) : einstein_coeff(di)
{
    char text_line[MAX_TEXT_LINE_WIDTH];
    int	i, i_max, parity, v, up, low;
    double coeff, energy, j, omega, hf;

    string file_name;
    ifstream input;

    if (verbosity)
        cout << "OH molecule radiative coefficients (including HF splitting) are being initializing..." << endl;

    file_name = path + "spectroscopy/radiative_oh_hf.txt";
    input.open(file_name.c_str(), ios_base::in);

    if (!input) {
        cout << "Error in " << SOURCE_NAME << ": can't open " << file_name << endl;
        exit(1);
    }
    input.getline(text_line, MAX_TEXT_LINE_WIDTH);
    input.getline(text_line, MAX_TEXT_LINE_WIDTH);
    input >> i_max;

    for (i = 0; i < i_max; i++)
    {
        input >> v >> j >> omega >> parity >> hf;
        up = di->get_nb(parity, v, j, omega, hf);

        input >> v >> j >> omega >> parity >> hf;
        low = di->get_nb(parity, v, j, omega, hf);

        input >> coeff >> energy;
        if (low != -1 && up != -1)
        {
            arr[up][low] = coeff;
            arr[low][up] = di->lev_array[up].g * coeff / ((double)di->lev_array[low].g);
        }
    }
    input.close();
    if (verbosity)
        cout << "  data are read from file " << file_name << endl;
}

nh3_einstein_coeff::nh3_einstein_coeff(const string &path, const energy_diagram *di, int verbosity) : einstein_coeff(di)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int	i, i_max, syminv, v, up, low;
	double coeff, energy, j, k;
	
	string file_name;
	ifstream input;
	
	if (verbosity) 
		cout << "NH3 molecule radiative coefficients are being initializing..." << endl;
	
	file_name = path + "spectroscopy/radiative_nh3.txt";
	input.open(file_name.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << file_name << endl;
		exit(1);
	}	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> i_max;
	
	for (i = 0; i < i_max; i++)
	{
		input >> v >> j >> k >> syminv;
		up = di->get_nb(syminv, v, j, k);
			
		input >> v >> j >> k >> syminv;
		low = di->get_nb(syminv, v, j, k);
			
		input >> coeff >> energy;
		if (low != -1 && up != -1)
		{
			arr[up][low] = coeff;
			arr[low][up] = di->lev_array[up].g *coeff / ((double) di->lev_array[low].g);
		}
	}
	input.close();
	if (verbosity) 
		cout << "  data are read from file " << file_name << endl;
}

h2co_einstein_coeff::h2co_einstein_coeff(const std::string& path, const energy_diagram *di, int verbosity) : einstein_coeff(di)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int	i, i_max, v, j, ka, kc, up, low;
	double coeff, energy;

	string file_name;
	ifstream input;

	if (verbosity)
		cout << "H2CO molecule radiative coefficients are being initializing..." << endl;

	file_name = path + "spectroscopy/radiative_h2co.txt";
	input.open(file_name.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << file_name << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> i_max;

	for (i = 0; i < i_max; i++)
	{
		input >> v >> j >> ka >> kc;
		up = di->get_nb(v, j, ka - kc);

		input >> v >> j >> ka >> kc;
		low = di->get_nb(v, j, ka - kc);

		input >> coeff >> energy;
		if (low != -1 && up != -1)
		{
			arr[up][low] = coeff;
			arr[low][up] = di->lev_array[up].g * coeff / ((double)di->lev_array[low].g);
		}
	}
	input.close();
	if (verbosity)
		cout << "  data are read from file " << file_name << endl;
}


//
// Transition class
//

transition::transition(energy_level lowl, energy_level upl)
: low_lev(lowl), up_lev(upl)
{
	energy = up_lev.energy - low_lev.energy;
	// the frequency in Hz, energy in cm-1, speed in cm/s
	freq = SPEED_OF_LIGHT *energy; 
}

//
// Functions
//

void boltzmann_populations(double gas_temp, double *arr, const energy_diagram *diagram)
{
	int i, nb_mol_lev;
	double stat_sum(0.);
	
	nb_mol_lev = diagram->nb_lev;
	for (i = 0; i < nb_mol_lev; i++)
	{
		arr[i] = diagram->lev_array[i].g *exp(-diagram->lev_array[i].energy *CM_INVERSE_TO_KELVINS/gas_temp);
		stat_sum += arr[i];
	}
	for (i = 0; i < nb_mol_lev; i++) {
		arr[i] /= stat_sum;
	}
}

// Calculation results for CH3OH: the depth is equal to 5.e+15 at 100 K and thermal velocity dispersion 1.e+5 cm/s;
double calc_depth(const energy_diagram *diagram, const einstein_coeff *einst_c, double gas_temp, double vel_disp)
{
	int i, j, nb_lev;
	double a, en, d, d_min;
	double *pop;

	nb_lev = einst_c->nb_lev;
	pop = new double [nb_lev];

	// Boltzmann populations of the molecule are considered;
	boltzmann_populations(gas_temp, pop, diagram);
	
	d_min = 1.e+99;
	for (i = 1; i < nb_lev; i++) {
		for (j = 0; j < i; j++)
		{
			if ((a = einst_c->arr[i][j]) > 0.) 
			{
				// the molecule concentration is taken to be 1 cm^{-3}, all values are in Gauss system units;
				en = diagram->lev_array[i].energy - diagram->lev_array[j].energy;
				d = en*en*en *EIGHT_PI *SQRT_PI *vel_disp 
					/(a*(diagram->lev_array[i].g *pop[j]/ ((double) diagram->lev_array[j].g) - pop[i])); 
				
				if (d < d_min) 
					d_min = d;
			}
		}
	}
	delete [] pop;
	return d_min; // the answer is in cm;
}

void save_populations(const string &output_path, const energy_diagram *diagram, double *arr, int nb_cloud_lay, int nb_mol_lev, 
					  bool normalized, string id)
{
	int i,j;
	double w;
	
	string file_name;
	ofstream text_file;

	file_name = output_path + diagram->mol.name;
	file_name += "_populations";
	file_name += id;
	file_name += ".txt";
	text_file.open(file_name.c_str());

	if (!text_file.is_open()) 
		cout << "Error in " << SOURCE_NAME << ": can't open file to save population values;" << endl;
	else 
	{
		text_file.setf(ios::scientific);
		text_file.precision(6);
		text_file << left << setw(8) << nb_mol_lev << setw(8) << nb_cloud_lay << endl;
		
		text_file << setw(6) << " ";
		for (j = 0; j < nb_cloud_lay; j++) {
			text_file << left << setw(12) << j;
		}
		text_file << endl;
		
		for (i = 0; i < nb_mol_lev; i++)
		{
			if (normalized) 
				w = 1./ ((double) diagram->lev_array[i].g);
			else w = 1.;

			text_file << left << setw(6) << i;
			for (j = 0; j < nb_cloud_lay; j++) {
				text_file << left << setw(12) << arr[j*nb_mol_lev+i] *w;
			}
			text_file << endl;
		}
	}
	text_file.close();
}

void read_populations(const string &file_name, double *arr, int nb_cloud_lay, int nb_mol_lev)
{
	int i, j, k;
	ifstream text_file;

	text_file.open(file_name.c_str());

	if (!text_file.is_open()) 
		cout << "Error in " << SOURCE_NAME << ": can't open file " << file_name << endl;
	else 
	{	
		text_file >> i >> j;
		if (i != nb_mol_lev || j != nb_cloud_lay) 
			cout << "Error in " << SOURCE_NAME << ": the numbers of columns and rows are not correct in " << file_name << endl;
		else
		{
			for (j = 0; j < nb_cloud_lay; j++) 
				text_file >> k;

			for (i = 0; i < nb_mol_lev; i++) 
			{
				text_file >> k;
				for (j = 0; j < nb_cloud_lay; j++) {
					text_file >> arr[j*nb_mol_lev+i];
				}
			}
		}
	}
	text_file.close();
}
