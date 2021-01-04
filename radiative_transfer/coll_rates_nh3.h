#pragma once
#include "spectroscopy.h"
#include "coll_rates.h"


// The ortho and para levels of NH3 do not inter-convert in inelastic collisions;
// data on NH3-He collisions, BASECOL data, 5 < T < 300 K;
// Machin & Roueff, J. Phys. B 38, p. 1519, 2005; ortho-NH3 j <= 7 (22 levels), para-NH3 - j <= 4 (16 levels);
class nh3_he_coll_data : public collision_data
{
public:
	nh3_he_coll_data(const std::string path, const energy_diagram *, int verbosity=1);
};

// data on NH3-H2 collisions; 10 < T < 200 K;
// Bouhafs et al., MNRAS 470, p.2204, 2017; ortho-NH3 - 17 levels, para-NH3 - 34 levels;
// p-H2 is thermalized over j=0,2
class nh3_ph2_coll_data : public collision_data
{
public:
	nh3_ph2_coll_data(const std::string path, const energy_diagram *, int verbosity=1);
};

// o-H2 is j=1;
class nh3_oh2_coll_data : public collision_data
{
public:
	nh3_oh2_coll_data(const std::string path, const energy_diagram *, int verbosity=1);
};

// data on NH3-H collisions, 10 < T < 200 K;
// Bouhafs et al., MNRAS 470, p.2204, 2017; ortho-NH3 - 17 levels, para-NH3 - 34 levels;
class nh3_h_coll_data : public collision_data
{
public:
	nh3_h_coll_data(const std::string path, const energy_diagram *, int verbosity=1);
};

// Methods of this class calculate the total collisional rates for transitions between levels of NH3;
class nh3_collisions : public collisional_transitions
{
public:
	void set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, double h_conc,
		double el_conc, double *&concentration, int *&indices) const;

	void get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, double &up_rate,
		double temp_neutrals, const double *concentration, const int *indices) const;

	nh3_collisions(const std::string &, const energy_diagram *, int verbosity =1);
};
