#pragma once
#include "spectroscopy.h"
#include "coll_rates.h"

// ortho- and para-H2CO do not inter-convert in inelastic collisions;
// data taken from the database LAMBDA, 
// Wiesenfeld, Faure, MNRAS 432, 2573 (2013); 40 levels of ortho-H2CO, 41 levels of para-H2CO; 10 <= T <= 300 K;
class h2co_h2_coll_data : public collision_data
{
public:
	h2co_h2_coll_data(const std::string path, const energy_diagram *, bool is_h2_ortho, int verbosity = 1);
};

// data taken from the database BASECOL,
// Green, Astroph. J. Suppl. Ser. 76, 979 (1991); 40 levels of ortho-H2CO, 41 levels of para-H2CO; 10 <= T <= 300 K;
class h2co_he_coll_data : public collision_data
{
public:
	h2co_he_coll_data(const std::string path, const energy_diagram *, int verbosity = 1);
};


// Methods of this class calculate the total collisional rates for transitions between levels of H2CO;
class h2co_collisions : public collisional_transitions
{
public:
	void set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, double h_conc,
		double el_conc, double*& concentration, int*& indices) const;

	void get_rate_neutrals(const energy_level& up_lev, const energy_level& low_lev, double& down_rate, double& up_rate,
		double temp_neutrals, const double* concentration, const int* indices) const;

	h2co_collisions(const std::string&, const energy_diagram*, int verbosity = 1);
};
