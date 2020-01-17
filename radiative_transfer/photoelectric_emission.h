#pragma once

#include <string>
#include "radiation_field.h"
// some arbitrary small value:
#define MIN_PHOTOEMISSION_RATE 1.e-99

// Calculation of photoelectron emission for all grain sizes;
void calc_grain_photoelectron_rates(const std::string &input_path, double cr_uv_flux);

// Methods of this class calculate the photoelectric rate for dust particles with given radius and charge;
// description of the calculations see in Weingartner & Draine, ApJ Suppl S 134, p.263 (2001);
class photoel_emission
{
protected:
	int index_r, nb_ph_en, nb_ph_en2, nb_gr_rad, verbosity;
	double radius, length_el, work_f, t; // auxulary variables, store important parameters;
	
	// ph_en_arr[] - abs_coeff[][], ph_en_arr2[] - abs_length[], energy in cm-1
	double *ph_en_arr, *ph_en_arr2, *grain_rad_arr, *abs_length;
	double **abs_coeff; // first variable - photon energy, second - grain radius;
	const radiation_field *rfield;

	// radius and length in cm, energy in cm-1, cross section in cm2
	virtual void set_radius(double r);	
	virtual double get_abs_cs(double en) const;
	double get_abs_length(double en) const;

public:
	void set_rfield(const radiation_field *rf) { rfield = rf; }
	// in the calculations of zmax, the maximal photon energy is assumed to be 13.6 eV;
	virtual void get_limit_charges(double radius, int & zmin, int & zmax) const = 0;

	// calculation of photoelectric yield, energy in cm-1:
	virtual double get_yield(int z, double radius, double eph) const = 0;

	// integration over energies 0 - 13.6 eV;
	// calculation of the rate of electron emission by one grain (s-1) and average energy of electron (cm-1), radius in cm:
	virtual void get_phel_emiss(int z, double radius, double &rate, double &av_energy) = 0;
	
	photoel_emission(int verbosity);
	virtual ~photoel_emission();
};

class photoel_emission_carbon_dust : public photoel_emission
{
private:
	int nb_ph_en3, nb_pah_rad;
	double *ph_en_arr3, *pah_rad_arr;
	double **abs_coeff_p; // data for pah molecules
	
	void set_radius(double r);	
	double get_abs_cs(double en) const;

public:
	void get_limit_charges(double radius, int & zmin, int & zmax) const;
	double get_yield(int z, double radius, double eph) const;
	void get_phel_emiss(int z, double radius, double &rate, double &av_energy);
	
	photoel_emission_carbon_dust(const std::string &path, int verbosity = 1);
	~photoel_emission_carbon_dust();
};

// For silicate grains, the data on absorption coeffisient exist for radii 10 A <= r <= 1e+5 A = 10 um
class photoel_emission_silicate_dust : public photoel_emission
{
public:
	void get_limit_charges(double radius, int & zmin, int & zmax) const;
	double get_yield(int z, double radius, double eph) const;
	void get_phel_emiss(int z, double radius, double &rate, double &av_energy);
	
	photoel_emission_silicate_dust(const std::string &path, int verbosity = 1);
};
