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
	// lay_nb_hg is the layer nb with the highest gain, if gain < 0 (no amplification) for the entire cloud, lay_nb_hg = 0; 
	int nb_cloud_lay, lay_nb_hg, nb_aspect_ratio, nb_freq;  
    // population inversion = n_up/g_up - n_l/g_l [cm-3], gain > 0 for amplification [cm-1], exc_temp < 0 [K]
    // inv, gain, lum - average value of the inversion, gain, luminosity over entire cloud;
	// tau_max - the optical depth summed over regions with positive gain along the outflow at line centre (taking into account the frequency shift)
	// tau_eff - the effective optical depth (ignoring frequency shift),
	// tau_sat - the saturation optical depth at the given point,
	double inv, gain, lum, tau_max, tau_eff, tau_sat;
	// loss rate is an average for the upper and lower levels,
    // pump efficiency = inversion/(n_up/g_up + n_l/g_l) - is dimensionless, loss rate in [s-1], luminosity in [cm-3 s-1] 
    double *inv_arr, *gain_arr, *lum_arr, *emiss_coeff_arr, *pump_rate_arr, *pump_eff_arr, *loss_rate_arr, *exc_temp_arr, 
		*tau_vs_aspect_ratio, *tau_vs_frequency;

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
	double  velocity_shift;
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
    // gain and tau are calculated, the nb of the layer with the highest gain is saved,
	void calc_gain(transition_data &, double *level_pop);

	// calculation of optical depth (taking into account the frequency shift) as a function of aspect ratio a = 1/cos theta,
	// theta - angle between line of sight and z,
	void calc_line_profile(transition_data&, double* level_pop);

	// the excitation temperature of the levels is calculated,
	void calc_exc_temp(transition_data&, double* level_pop);
	
    // calculates parameter tau_sat,
    // the parameter loss_rate and excitation temperature must be calculated (the functions lim_luminosity_lvg(), calc_gain(), calc_exc_temp() must be called before),
	// beaming factor is dOmega/4pi, for isotropic is equal to 1.
    void calc_saturation_depth(double beaming_factor);
    
	// the function deletes the old data on inverted transitions, must be called first,    
    // the inversion, gain, optical depth and excitation temperature are calculated,
    // the transition is added if:
    // 1. the population inversion in one layer > population error; 2. optical depth is > min_optical_depth,
	void find(double *level_pop, double rel_error);

    // must be called after the previous function; the inversion, gain, optical depth and excitation temperature are calculated,
    void add(std::list<transition>& trans_list, double* level_pop); 

	const transition_data* get(const transition *) const;
	
	void save_short_data(int process_nb, const std::string & fname) const;
	void save_full_data(const std::string & fname) const;
    void save_transition(const std::string& fname, const transition& trans) const;
	
	// the dependence of optical depth on aspect ratio,
	void save_optical_depth_1(const std::string& fname) const; 
	// the dependence of optical depth along the shock outflow on velocity (frequency)
	void save_optical_depth_2(const std::string& fname) const;
	
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
