#pragma once

#include <list>
#include <cfloat>
#include <cmath>
#include "constants.h"
#include "spectroscopy.h"
#include "cloud_data.h"
#include "coll_rates.h"

class transition_data
{
public:
	int nb_cloud_lay;
    // population inversion = n_up/g_up - n_l/g_l [cm-3], gain > 0 for amplification [cm-1], exc_temp < 0 [K]
    // average value of the gain over entire cloud, optical depth summed over regions with positive gain, 
	double inv, gain, lum, tau, tau_sat;  
    // pump efficiency = inversion/(n_up/g_up + n_l/g_l) - is dimensionless, loss rate in [s-1], luminosity in [cm-3 s-1] 
    double *inv_arr, *gain_arr, *lum_arr, *pump_eff_arr, *loss_rate_arr, *exc_temp_arr;

	const transition *trans;
	
	// The assignment and relation operators are needed to sort by energy
	transition_data& operator = (const transition_data &obj);
	bool operator < (const transition_data &obj) { return (*trans < *obj.trans); };
	bool operator > (const transition_data &obj) { return (*trans > *obj.trans); };
	bool operator != (const transition_data &obj) { return *trans != *obj.trans; };
	bool operator == (const transition_data &obj) { return *trans == *obj.trans; };
	
	transition_data(int nb_cloud_lay, const energy_level &low, const energy_level &up);
	transition_data(const transition_data &);
	virtual ~transition_data();
};

class transition_data_container
{
protected:
	int		t_nb;
    double  min_optical_depth;  // it is necessary for the function find(), by default is very low
	double	*tau_arr, *sm_arr;

	const energy_diagram	*diagram;
	const einstein_coeff	*einst_coeff;
	const cloud_data		*cloud;
	
public:
	int nb_mol_lev, nb_cloud_lay;
	std::list<transition_data> data;
	
    void set_min_optical_depth(double depth) { min_optical_depth = depth; }
	void calc_inv(transition_data &, double *level_pop);
    
    // check water molecule name in this routine, the calculated gain takes into account the absorption by dust,
    // gain and tau are calculated,
	void calc_gain(transition_data &, double *level_pop);

	// the excitation temperature of the levels is calculated;
	void calc_exc_temp(transition_data&, double* level_pop);
	
    // calculates parameter tau_sat,
    // the parameter loss_rate must be calculated (the function lim_luminosity_lvg must be called before),
    void calc_saturation_depth(double beaming_factor);
    
	// the function deletes the old data on inverted transitions, must be called first,    
    // the inversion, gain and excitation temperature are calculated,
    // the transition is added if:
    // 1. the population inversion in one layer > population error; 2. optical depth is > min_optical_depth,
	void find(double *level_pop, double rel_error);

    // must be called after the previous function; the inversion, gain and excitation temperature are calculated,
    void add(std::list<transition>& trans_list, double* level_pop); 

	const transition_data* get(const transition *) const;
	
	void save_short_data(int process_nb, const std::string & fname) const;
	void save_full_data(const std::string & fname) const;
    void save_transition(const std::string& fname, const transition& trans) const;
	
	virtual ~transition_data_container();
	transition_data_container(const cloud_data *, const energy_diagram *, const einstein_coeff *);
};

struct saturated_maser_func
{
	double tau; // tau is the optical depth in the line center;
	double operator() (double u) {
		return ONEDIVBY_SQRT_PI * exp(tau*exp(-u*u) - u*u);
	}
};
