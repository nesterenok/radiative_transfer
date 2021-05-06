
#pragma once
#include <string>
#include <list>
#include <vector>
#include "spectroscopy.h"

// some arbitrary small value of the rate is used, the rates lower than this constant are set to zero:
#define MIN_COLLISION_RATE 1.e-99

// The class with the data on collisional coefficients;
class collision_data
{
public:
	int		imax, jmax;	// the dimensions of coeff[][] array; jmax is the dimension of tgrid[] array;
	int		nb_lev;		// number of molecule levels for which the coefficients are available, imax = nb_lev*(nb_lev-1)/2;
	double	*tgrid;		// the array with temperature grid points;
	// the data in array are arranged 1->0, 2->0, 2->1, 3->0,..
	// the dimension of coeff_deriv[][] must be imax*jmax in the case of cubic spline implementation;
	double	**coeff, **coeff_deriv; 

// return the maximal temperature in the data set;
	double get_max_temp() const { return tgrid[jmax-1]; }

// these functions find the nb of the lower boundary of the interval l, such that tgrid[l] < temp < tgrid[l+1];
	int locate(double temp) const; 
    int hunt_index(double temp, int old_index) const; //

// the function calculates the array of rate derivatives;
	virtual void calc_coeff_deriv();

// in the functions below the linear approximation of the data is employed, first_lev > sec_lev;
// there is no check that the level numbers and temperature are in the ranges of the arrays coeff[][] and tgrid[], 
// this check must be done in the function, that calls this method;
	virtual double get_rate(int first_lev, int sec_lev, double temp) const;
// for the fast work the index for the tgrid[] array must given here, tgrid[lo] < temp < tgrid[lo+1]: 
	virtual double get_rate(int first_lev, int sec_lev, int lo, double temp) const;
	
	collision_data();
	virtual ~collision_data(); 
};

// may be it is better to init own array of spline coefficients...
class collision_data_cub_spline : public collision_data
{
public:
	// the function calculates the array of second order derivatives of the data;
	void calc_coeff_deriv();
	
	// the index lo must be < jmax-1 here;
	// for low rate values (low temperatures), the interpolation may give negative rates;
	double get_rate(int first_lev, int sec_lev, int lo, double temp) const;
};

// The class, that calculates the rates of collisional transitions;
class collisional_transitions
{
protected:
	int	nb_lev, nb1, nb2, nb3; // nb3 is the total nb of collision data;
	// Maximal temperatures of the collisional data:
	double *max_temp;
	std::vector<collision_data *> coll_data;

public:
	// arrays with concentration of collisional partners and indices are initialized in the routine:
	virtual void set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, double h_conc, 
        double el_conc, double *&concentration, int *&indices) const;

    virtual void set_ion_param(double temp_neutrals, double temp_ions, double hp_conc, double h3p_conc, 
        double *&concentration, int *&indices) const;

	// the function calculates the collision rates due to neutral species:
	virtual void get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, double &up_rate, 
		double temp_neutrals, const double *concentration, const int *indices) const;
	
	// the function calculates the collision rates due to electrons:
	virtual void get_rate_electrons(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, double &up_rate, 
		double temp_el, const double *concentration, const int *indices) const;
    
    // the function calculates the collision rates due to ions:
    virtual void get_rate_ions(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, double &up_rate,
        double temp_neutrals, double temp_ions, const double *concentration, const int *indices) const;


	double get_rate_neutrals(const energy_level &init_lev, const energy_level &fin_lev, 
		double temp_neutrals, const double *concentration, const int *indices) const;
	
	double get_rate_electrons(const energy_level &init_lev, const energy_level &fin_lev, 
		double temp_el, const double *concentration, const int *indices) const;

    double get_rate_ions(const energy_level &init_lev, const energy_level &fin_lev,
        double temp_neutrals, double temp_ions, const double *concentration, const int *indices) const;

    // initial level il, final level fl, il > fl;
    void check_spline(int il, int fl, std::string path, std::string name) const;

	collisional_transitions();
	virtual ~collisional_transitions(); // all collisional data are deleted here
};

// Dissociation data
class dissociation_data
{
protected:
	int jmax; // the dimension of tgrid[] array;
	int imax; // number of molecule levels for which the coefficients are available,
	double *tgrid; // the array with temperature grid points;
	double **coeff, **coeff_deriv;

public:
	// i - number of the level, temp - gas temperature, linear interpolation is used;
    // there is linear extrapolation of the rate at temperatures higher than the maximal value for which data exist,
	virtual double get_rate(int i, double temp) const;
	
	// the function calculates the array of rate derivatives;
	virtual void calc_coeff_deriv();
	
	dissociation_data();
	virtual ~dissociation_data();
}; 


// The function calculates the populations for optically thin media;
// the arrays with concentrations of collisional partners and indices have to be initialized before this function call; 
void opt_thin_pop(double *popul, const energy_diagram *, const einstein_coeff *, const collisional_transitions *, 
	double temp_neutrals, double temp_el, const double *concentration, const int *indices);

// Calculation of the heating rate of the gas components; by definition cooling rate is < 0., heating rate > 0.;
double heating_of_neutral_gas(double *arr, const energy_diagram *diagram, const collisional_transitions *coll_trans,
	double temp_neutrals, const double *concentration, const int *indices);

double heating_of_electron_gas(double *arr, const energy_diagram *diagram, const collisional_transitions *coll_trans, 
	double temp_el, const double *concentration, const int *indices);
