
/*	Modifications:
	11.05.2016. The class collision_data_cub_spline was added. The cubic spline interpolation of rate coefficients were added. 
	12.05.2016. Index l_h was added to the class collisional_transitions;
	27.05.2016. The function opt_thin_pop() was added. It is the analogue of boundary_layer_populations(), and, may be,
		has to replace it;
	22.06.2016. The destruction of arrays were added to the destructor of base class collision_data(); 
		Note: the deletion of the data must be removed from derived classes!
		Electron temperature and index were added to the parameter list of the class collisional_transitions;
	27.06.2016. Significant changes to the variable list in the class collisional_transitions. The scheme for the calculation
		of rate coefficients was modified.
	04.03.2017. Check for errors. The error was found in the formula for cubic spline interpolation.
	08.03.2017. Check for errors. 
*/

#include <stdlib.h>
#include <memory.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <sstream>

#include "constants.h"
#include "linear_algebra.h"
#include "coll_rates.h"

using namespace std;
const double one_sixth_const = 1. / 6.;

//
// The classes that contain rate coefficient data
//

collision_data::collision_data() 
: imax(0), jmax(0), nb_lev(0), tgrid(0), coeff(0), coeff_deriv(0)
{;}

collision_data::~collision_data()
{
	delete [] tgrid;
	if (coeff != 0) 
		free_2d_array(coeff);
	
	if (coeff_deriv != 0) 
		free_2d_array(coeff_deriv);
}

// The energy of the first level must be higher, first_lev > sec_lev
double collision_data::get_rate(int first_lev, int sec_lev, double temp) const {
	return get_rate(first_lev, sec_lev, locate(temp), temp);
}

double collision_data::get_rate(int first_lev, int sec_lev, int lo, double temp) const
{
	// the data in array are arranged 1->0, 2->0, 2->1, 3->0,..
	int i = first_lev*(first_lev -1)/2 + sec_lev;
	return coeff[i][lo] + coeff_deriv[i][lo] *(temp - tgrid[lo]);
}

void collision_data::calc_coeff_deriv()
{
	int i, j;
	for (i = 0; i < imax; i++) {
		for (j = 0; j < jmax-1; j++) {
			coeff_deriv[i][j] = (coeff[i][j+1] - coeff[i][j])/(tgrid[j+1] - tgrid[j]);
		}
	}
}

int collision_data::locate(double temp) const
{
	int j, l = 0, r = jmax-1; 
	while (r-l > 1)
	{
		j = l + ((r-l) >> 1);
		if (tgrid[j] < temp) 
			l = j;
		else r = j;
	}
	return l;
}

int collision_data::hunt_index(double temp, int old_index) const
{
    return 0;
}

void collision_data_cub_spline::calc_coeff_deriv()
{
	int i, j;
	double p, sig;
	double *u = new double [jmax-1];
	
	for (i = 0; i < imax; i++) 
	{
	// the lower boundary condition is set to be "natural"
		coeff_deriv[i][0] = u[0] = 0.;
	
		// this is the decomposition loop of the tridiagonal algorithm. coeff_deriv[][] and u are used for temporary
		// storage of the decomposed factors.
		for (j = 1; j < jmax-1; j++) 
		{ 
			sig = (tgrid[j] - tgrid[j-1])/(tgrid[j+1] - tgrid[j-1]);
			p = sig *coeff_deriv[i][j-1] + 2.;

			coeff_deriv[i][j] = (sig - 1.)/p;
			u[j] = (coeff[i][j+1] - coeff[i][j])/(tgrid[j+1] - tgrid[j]) - (coeff[i][j] - coeff[i][j-1])/(tgrid[j] - tgrid[j-1]);
			u[j] = (6.*u[j]/(tgrid[j+1] - tgrid[j-1]) - sig*u[j-1])/p;
		}
		// the upper boundary condition is set to be "natural"
		coeff_deriv[i][jmax-1] = 0.;
	
		for (j = jmax-2; j >= 0; j--) { // this is the backsubstitution loop of the tridiagonal algorithm.
			coeff_deriv[i][j] = coeff_deriv[i][j]*coeff_deriv[i][j+1] + u[j]; // 
		}
	}
	delete [] u;
}

double collision_data_cub_spline::get_rate(int first_lev, int sec_lev, int lo, double temp) const
{
	int i;
	double a, b, h;
	
	i = ((first_lev*(first_lev -1)) >> 1) + sec_lev;
	h = tgrid[lo+1] - tgrid[lo];
	a = (tgrid[lo+1] - temp)/h;
	b = 1. - a; // = (temp - tgrid[lo])/h;
	
	return fabs(a*coeff[i][lo] + b*coeff[i][lo+1] + one_sixth_const*((a*a*a - a)*coeff_deriv[i][lo] + (b*b*b - b)*coeff_deriv[i][lo+1])*(h*h));
}

//
// The classes that calculate collisional rates
//

collisional_transitions::collisional_transitions()
: nb_lev(0), nb1(0), nb2(0), nb3(0), max_temp(0)
{;}

collisional_transitions::~collisional_transitions()
{
	for (int i = 0; i < nb3; i++) {
		delete coll_data[i];
	}
	delete [] max_temp;
}

// The function returns arrays with indices and concentrations of collisional partners in order to speed up the calculations;
// Previous values of the pointers are ignored.
void collisional_transitions::set_gas_param(double temp_neutrals, double temp_el, double hec, double ph2c, double oh2c, double hc, 
    double ec, double *&concentration, int *&indices) const
{
	int i;
	if (indices != 0)
		delete [] indices;
	
	if (concentration != 0)
		delete [] concentration;

	indices = new int [nb3];
	concentration = new double [nb3];

	// neutral collisional partners:
	for (i = 0; i < nb1; i++) {
		indices[i] = coll_data[i]->locate(temp_neutrals);
	}
	// collisions with electrons:
	for (i = nb1; i < nb2; i++) {
		concentration[i] = ec;
		indices[i] = coll_data[i]->locate(temp_el);
	}
}

void collisional_transitions::set_ion_param(double temp_neutrals, double temp_ions, double hp_conc, double h3p_conc,
    double *&concentration, int *&indices) const
{;}

// The energy of the first level is higher, up_lev.nb > low_lev.nb;
void collisional_transitions::get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate,
    double &up_rate, double temp_neutrals, const double *concentration, const int *indices) const
{
    down_rate = 0.;
    for (int i = 0; i < nb1; i++)
    {
        if (up_lev.nb < coll_data[i]->nb_lev) {
            down_rate += coll_data[i]->get_rate(up_lev.nb, low_lev.nb, indices[i], (temp_neutrals < max_temp[i]) ? temp_neutrals : max_temp[i])
                *concentration[i];
        }
    }

    if (down_rate > MIN_COLLISION_RATE) {
        up_rate = down_rate * exp((low_lev.energy - up_lev.energy)*CM_INVERSE_TO_KELVINS / temp_neutrals) *up_lev.g / ((double)low_lev.g);
    }
    else down_rate = up_rate = 0.;
}

void collisional_transitions::get_rate_electrons(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, 
	double &up_rate, double temp_el, const double *concentration, const int *indices) const
{
	down_rate = 0.;
	// there may be more than one data set for collisions with electrons, the data must be ordered by their significance;
	for (int i = nb1; i < nb2; i++) 
	{
		if (up_lev.nb < coll_data[i]->nb_lev) {
			down_rate = coll_data[i]->get_rate(up_lev.nb, low_lev.nb, indices[i], (temp_el < max_temp[i]) ? temp_el : max_temp[i])
				*concentration[i];
			break;
		}
	}

	if (down_rate > MIN_COLLISION_RATE) {
		up_rate = down_rate *exp((low_lev.energy - up_lev.energy)*CM_INVERSE_TO_KELVINS/temp_el) *up_lev.g/((double) low_lev.g);
	}
	else down_rate = up_rate = 0.;
}

void collisional_transitions::get_rate_ions(const energy_level &up_lev, const energy_level &low_lev, double &down_rate,
    double &up_rate, double temp_neutrals, double temp_ions, const double *concentration, const int *indices) const
{
   down_rate = up_rate = 0.;
}

double collisional_transitions::get_rate_neutrals(const energy_level &init_lev, const energy_level &fin_lev, 
	double temp_neutrals, const double *concentration, const int *indices) const
{
	double down, up;
	if (init_lev.nb > fin_lev.nb)
	{
		get_rate_neutrals(init_lev, fin_lev, down, up, temp_neutrals, concentration, indices);
		return down;
	}
	else if (init_lev.nb < fin_lev.nb)
	{
		get_rate_neutrals(fin_lev, init_lev, down, up, temp_neutrals, concentration, indices);
		return up;
	}
	return 0.;
}

double collisional_transitions::get_rate_electrons(const energy_level &init_lev, const energy_level &fin_lev,
	double temp_el, const double *concentration, const int *indices) const
{
	double down, up;
	if (init_lev.nb > fin_lev.nb)
	{
		get_rate_electrons(init_lev, fin_lev, down, up, temp_el, concentration, indices);
		return down;
	}
	else if (init_lev.nb < fin_lev.nb)
	{
		get_rate_electrons(fin_lev, init_lev, down, up, temp_el, concentration, indices);
		return up;
	}
	return 0.;
}

double collisional_transitions::get_rate_ions(const energy_level &init_lev, const energy_level &fin_lev,
    double temp_neutrals, double temp_ions, const double *concentration, const int *indices) const
{
    double down, up;
    if (init_lev.nb > fin_lev.nb)
    {
        get_rate_ions(init_lev, fin_lev, down, up, temp_neutrals, temp_ions, concentration, indices);
        return down;
    }
    else if (init_lev.nb < fin_lev.nb)
    {
        get_rate_ions(fin_lev, init_lev, down, up, temp_neutrals, temp_ions, concentration, indices);
        return up;
    }
    return 0.;
}

void collisional_transitions::check_spline(int il, int fl, std::string path, std::string name) const
{
    int nb = 5;
    double a, rate, t;
    stringstream ss;
    ofstream output;

    ss.clear();
    ss << path << name << "_coll_coeff_" << il << "_" << fl << ".txt";

    output.open(ss.str().c_str(), ios_base::out);
    output << scientific;
    output.precision(3);

    t = 10.;
    output << left << setw(12) << "! temp(K)";
    for (nb = 0; nb < (int)coll_data.size(); nb++) {
        output << left << setw(12) << nb;
    }
    output << endl;

    for (t = 10; t < 10000.; t *= 1.1) {
        output << left << setw(12) << t;

        for (nb = 0; nb < (int)coll_data.size(); nb++) {
            if (il < coll_data[nb]->nb_lev) {
                a = coll_data[nb]->get_max_temp();
                rate = (t < a) ? coll_data[nb]->get_rate(il, fl, t) : coll_data[nb]->get_rate(il, fl, a);
                output << left << setw(12) << rate;
            }
        }
        output << endl;
    }
    output.close();
}

//
// The classes that contain data on dissociation rates 
//

dissociation_data::dissociation_data() 
	: imax(0), jmax(0), tgrid(0), coeff(0), coeff_deriv(0)
{;}

dissociation_data::~dissociation_data() 
{
	delete [] tgrid;
	if (coeff != 0) 
		free_2d_array(coeff);
	
	if (coeff_deriv != 0) 
		free_2d_array(coeff_deriv);
}

void dissociation_data::calc_coeff_deriv()
{
	int i, j;
	for (i = 0; i < imax; i++) {
		for (j = 0; j < jmax-1; j++) {
			coeff_deriv[i][j] = (coeff[i][j+1] - coeff[i][j])/(tgrid[j+1] - tgrid[j]);
		}
	}
}

// the linear extrapolation is used here, check for maximal speed is not implemented for dissociation,
double dissociation_data::get_rate(int i, double temp) const 
{
//  if (temp > tgrid[jmax-1])
//      return coeff[i][jmax-1];

	int j, l = 0, r = jmax-1; 
	while (r-l > 1)
	{
		j = l + ((r-l) >> 1);
		if (tgrid[j] < temp) 
			l = j;
		else r = j;
	}
	return coeff[i][l] + coeff_deriv[i][l] *(temp - tgrid[l]);
}

//
// Functions
//

void opt_thin_pop(double *arr, const energy_diagram *diagram, const einstein_coeff *e_coeff, const collisional_transitions *coll_trans,
	double temp_neutrals, double temp_el, const double *concentration, const int *indices)
{
	int i, j, nb_lev;
	double down1, down2, up1, up2;
	double **matrix;

	nb_lev = diagram->nb_lev;
	matrix = alloc_2d_array<double>(nb_lev, nb_lev);
	
	memset(*matrix, 0, nb_lev*nb_lev*sizeof(double));
	memset(arr, 0, nb_lev*sizeof(double));
	
	for (i = 1; i < nb_lev; i++) {
		for (j = 0; j < i; j++)
		{
			coll_trans->get_rate_neutrals(diagram->lev_array[i], diagram->lev_array[j], down1, up1,
				temp_neutrals, concentration, indices);

			coll_trans->get_rate_electrons(diagram->lev_array[i], diagram->lev_array[j], down2, up2,
				temp_el, concentration, indices);

			matrix[j][i] = e_coeff->arr[i][j] + down1 + down2;	// i->j
			matrix[i][i] -= e_coeff->arr[i][j] + down1 + down2;
		
			matrix[i][j] = up1 + up2;	// j->i
			matrix[j][j] -= up1 + up2;
		}
	}

	for (i = 0; i < nb_lev; i++) {
		matrix[0][i] = 1.;
	}
	arr[0] = 1.;
	
	lu_matrix_solve(matrix, arr, nb_lev);
	free_2d_array<double>(matrix);
}

double heating_of_neutral_gas(double *arr, const energy_diagram *diagram, const collisional_transitions *coll_trans, 
	double temp_neutrals, const double *concentration, const int *indices)
{
	int i, j;
	double down, up, heating(0.);

	for (i = 1; i < diagram->nb_lev; i++) {
		for (j = 0; j < i; j++)
		{
			coll_trans->get_rate_neutrals(diagram->lev_array[i], diagram->lev_array[j], down, up, 
				temp_neutrals, concentration, indices);
			heating += (down*arr[i] - up*arr[j]) *(diagram->lev_array[i].energy - diagram->lev_array[j].energy);
		}
	}
	// if arr is relative populations, answer is in erg s-1, and must be multiplied by the concentration of the coolant;
	return heating*CM_INVERSE_TO_ERG; 
}

double heating_of_electron_gas(double *arr, const energy_diagram *diagram, const collisional_transitions *coll_trans,
	double temp_el, const double *concentration, const int *indices)
{
	int i, j;
	double down, up, heating(0.);
	
	for (i = 1; i < diagram->nb_lev; i++) {
		for (j = 0; j < i; j++)
		{
			coll_trans->get_rate_electrons(diagram->lev_array[i], diagram->lev_array[j], down, up,
				temp_el, concentration, indices);
			heating += (down*arr[i] - up*arr[j]) *(diagram->lev_array[i].energy - diagram->lev_array[j].energy);
		}
	}
	return heating*CM_INVERSE_TO_ERG;
}
