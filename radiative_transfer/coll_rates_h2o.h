#pragma once
#include <string>
#include "spectroscopy.h"
#include "coll_rates.h"

// The function rewrites data files with collisional rate coefficients for h2o-e and h2o-h2 collisions by Faure & Josselin (2008);
// in original data files, the levels are labelled by number, but must to be labelled by (v,j,k1-k2);
void reformat_h2o_coll_data(const std::string &, int spin);

// The rate coefficients for the collisional transitions between lowest 45 levels of ortho(para) h2o molecule in h2o-h2 collisions; 
// A. Faure et al., Astron. Astrophys. 472, 1029 (2007).
// Additional point at 0 K was added, the data were set to zero at this temperature (the same is true for other data);
// Note, for the cold cloud simulations the low temperature data are needed;
class h2o_oh2_coll_data : public collision_data
{
public:
	h2o_oh2_coll_data(const std::string &, const energy_diagram *, int verbosity=1);
};

// A. Faure et al., Astron. Astrophys. 472, 1029 (2007).
class h2o_ph2_coll_data : public collision_data
{
public:
	h2o_ph2_coll_data(const std::string &, const energy_diagram *, int verbosity=1);
};

// The rate coefficients for the rovibrational transitions between 411 levels of ortho-h2o molecule and
// 413 levels of para-h2o molecule in h2o-h2 collisions; A. Faure and E. Josselin, Astron. Astrophys. 492, 257 (2008).
class h2o_h2_coll_rovibr_data : public collision_data
{
public:
	h2o_h2_coll_rovibr_data(const std::string &, const energy_diagram *, int verbosity=1);
};

// The rate coefficients for the collisional transitions between lowest 45 levels of ortho(para) h2o molecule
// in h2o-he collisions; S. Green et al., Astrophys. J. Suppl. Ser. 85, 181 (1993).
class h2o_he_coll_data : public collision_data
{
public:
	h2o_he_coll_data(const std::string &, const energy_diagram *, int verbosity=1);
};

// Extrapolated data of Green et al. (1993) to higher rovibrational levels of h2o, 411 levels of ortho(para) h2o;
// Nesterenok A.V., Astronomy Letters 39, pp. 717-728 (2013);
class h2o_he_coll_rovibr_data : public collision_data
{
public:
	h2o_he_coll_rovibr_data(const std::string &, const energy_diagram *, bool is_scaled, int verbosity=1);
};

// Rate coefficients for the rovibrational transitions between 411 levels of ortho-h2o and 413 levels of para-h2o in h2o-e collisions; 
// A. Faure and E. Josselin, Astron. Astrophys. 492, 257 (2008).
class h2o_e_coll_rovibr_data : public collision_data
{
public:
	h2o_e_coll_rovibr_data(const std::string &, const energy_diagram *, int verbosity=1);
};

// Rate coefficients for transitions betwen 45 levels of ortho/para-H2O and H atom
// Daniel et al., MNRAS 446, p. 2312, 2015;
class h2o_h_coll_data : public collision_data
{
public:
	h2o_h_coll_data(const std::string &, const energy_diagram *, int verbosity=1);
};


// The methods of this class calculate the total collisional rates for transitions between levels of H2O;
class h2o_collisions : public collisional_transitions
{
public:
	void set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, double h_conc,
		double el_conc, double *&concentration, int *&indices) const;

	void get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, double &up_rate,
		double temp_neutrals, const double *concentration, const int *indices) const;

	h2o_collisions(const std::string &, const energy_diagram *, bool he_is_scaled = false, int verbosity =1);
};
