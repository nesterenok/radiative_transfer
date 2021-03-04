#pragma once
#include <string>
#include "spectroscopy.h"
#include "coll_rates.h"

// Nb of torsionally excited states, for which the collisional rate coefficients were calculated in Rabli, Flower (2010, 2011)
#define NB_VIBR_EXCIT_CH3OH_RABLI 2
#define USE_TEMPER_EXTRAP_CH3OH 0  // 0 - false, 1 - true;

// The class includes data on rotationally and torsionally inelastic scattering of CH3OH-He system;
// Rabli, Flower, MNRAS 403, p.2033, 2010; Rabli, Flower, MNRAS 411, p.2093, 2011;
class ch3oh_he_coll_data : public collision_data
{
public:	
	ch3oh_he_coll_data(const std::string &, const energy_diagram *, int verbosity =1);
};

// Includes data on rotationally inelastic scattering of CH3OH and para-H2 system (Rabli, Flower, MNRAS 406, p.95, 2010);
class ch3oh_ph2_coll_data : public collision_data
{
public:	
	ch3oh_ph2_coll_data(const std::string &, const energy_diagram *, int verbosity =1);
};

// Includes data on rotationally inelastic scattering of CH3OH and ortho-H2 system (Rabli, Flower, MNRAS 406, p.95, 2010);
class ch3oh_oh2_coll_data : public collision_data
{
public:	
	ch3oh_oh2_coll_data(const std::string &, const energy_diagram *, int verbosity =1);
};

// The methods of this class calculate the total collisional rates for transitions between levels of CH3OH;
class ch3oh_collisions : public collisional_transitions
{
public:
	void set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, double h_conc,
		double el_conc, double *&concentration, int *&indices) const;
	
	void get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, double &up_rate,
		double temp_neutrals, const double *concentration, const int *indices) const;

	ch3oh_collisions(const std::string &, const energy_diagram *, int verbosity =1);
};
