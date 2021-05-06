#pragma once
#include <string>
#include <vector>
#include "photoelectric_emission.h"

#define NB_OF_DUST_MODELS_WD 3

// Draine & Li, ApJ 551, p.807, 2001;
#define GRAPHITE_DEBYE_TEMP 950 // in K
#define SILICATE_DEBYE_TEMP 500 

// density of dust material (for large grains), graphite - 2.24 g/cm3, silicate - 3.5 g/cm3, in g/cm3,
// the models on dust optical properties by prof. B.T. Draine are used, 
// the dust material density should be matched with the values given in data files;
#define GRAPHITE_MATERIAL_DENSITY 2.24 
#define SILICATE_MATERIAL_DENSITY 3.5

// size of large dust grains and PAH molecules, in cm:
#define RADIUS_OF_LARGE_GRAINS 0.5e-5
#define RADIUS_OF_PAH_MOLECULES 4.e-8 // Draine, Li, ApJ 657, pp. 810-837 (2007); PAH radius must be < 20 A; 1 A = 1.e-8 cm;


// material of large dust grains, 1 - silicate, 0 - carbonaceous, 
// this parameter determines the data on photon absorption cross section, photoelectric yields;
#define MATERIAL_OF_LARGE_DUST_GRAINS 1

// File names for the input data and output must be given. For carbonaceous and silicate grains.
// This routine rewrites the data file written in the format by prof. B. Draine to the format used by this code;
void reformat_dust_files(const std::string &path, std::string name1, std::string name2);

// For PAH carbonaceous grains,
// first file with data on PAH molecules, second - on graphite grains, third - output file;
void reformat_dust_files_PAH(const std::string &path, std::string name1, std::string name2, std::string name3);

// Grain distribution is normalized on total H nuclei concentration, f(a) = (1/n_H)*dn/da;
class size_distribution
{
public:
	// true - the grain size distribution; false - grains have one size; 
	bool distribution;
	double norm, rad, rad_min, rad_max;  // radius in cm;
	
	// there is no check if radius lies in the given range in the operator():
	virtual double operator() (double radius) const { return 0.; }
	virtual double operator() () const { return norm; }

	// for mono-size distribution:
	size_distribution(double r);
	// it is assumed that rmin < rmax:
	size_distribution(double rmin, double rmax);
};

class size_distribution_mono : public size_distribution
{
public:
	// parameters: radius and concentration per H;
	size_distribution_mono(double r, double conc);
};

const std::string WD_model_names[NB_OF_DUST_MODELS_WD] = {
	"W&D2001_Rv3.1_bc6e-5_A", "W&D2001_Rv5.5_bc0_B", "W&D2001_Rv5.5_bc3e-5_B"
};

// The table given by Weingartner & Draine, ApJ 548, p.296 (2001);
// All grain abundances must be rescaled relative to size distribution by Weingartner & Draine (2001) - this is done in the code;
// Rv = 3.1 scaling factor 0.93; Rv = 4.0 factor 1.18; Rv = 5.5 factor 1.42 (Draine, ARA&A 41, p. 241, 2003, http://www.astro.princeton.edu/~draine/dust/dustmix.html);
// for silicate grains, absorption data in the file is given for grain radii > 10 A = 1.e-7 cm
// first two parameters: Rv and bc; last four - minimal and maximal radii of carbonaceous grains, and for silicate grains (in cm);
const double WD2001_model_data[NB_OF_DUST_MODELS_WD][15] = {
	{3.1, 6.e-5, -1.54, -0.165, 0.0107, 0.428, 9.99e-12, -2.21, 0.3,    0.164, 1.e-13,   3.5e-8, 1.e-4, 1.e-7, 3.5e-5}, // Rv = 3.1; bc = 6.e-5;
	{5.5, 0.,    -2.8,  0.0356, 0.0203, 3.43,  2.74e-12, -1.09, -0.37,  0.218, 1.17e-13, 3.5e-8, 6.e-4, 5.e-7, 4.e-5},	 // Rv = 5.5; case B, b_c = 0;
	{5.5, 3.e-5, -1.9, -0.0517, 0.012,  7.28,  2.86e-12, -1.13, -0.109, 0.211, 1.04e-13, 3.5e-8, 1.e-3, 5.e-7, 4.e-5}};  // Rv = 5.5; case B, b_c = 3.e-5

// The grain size distribution according to Weingartner & Draine, ApJ 548, p.296 (2001); Draine & Li, ApJ 657, p.810 (2007);
class size_distribution_Draine2001 : public size_distribution
{
protected:
	double at, ac, alpha, beta, scaling_factor;

public:
	virtual double operator() (double radius) const;
	size_distribution_Draine2001(int model_nb, double rmin, double rmax);
};

// The size distribution for carbonaceous grains (including lognormal terms);
class size_distribution_carbon_Draine2001 : public size_distribution_Draine2001
{
protected:
	// carbon_abundance is the abundance of carbon in small grains (C/H ratio), = 0.93*6.e-5 for default distribution;
	double a1, a2, b1, b2, sigma1, sigma2, carbon_abundance;
	
public:
	// grain radius is in cm; 
	double operator() (double radius) const;
	size_distribution_carbon_Draine2001(int model_nb, double rmin, double rmax);
};

// The size distribution of silicate grains;
class size_distribution_silicate_Draine2001 : public size_distribution_Draine2001
{
public:
	size_distribution_silicate_Draine2001(int model_nb, double rmin, double rmax);
};

// Auxulary functors used in integration;
struct sd_func0
{
public:
	size_distribution *sd;
	double operator() (double radius) const { return (*sd)(radius); }
};

struct sd_func1
{
public:
	size_distribution *sd;
	double operator() (double radius) const { return (*sd)(radius)*radius; }
};

struct sd_func2
{
public:
	size_distribution *sd;
	double operator() (double radius) const { return (*sd)(radius)*radius*radius; }
};

struct sd_func3
{
public:
	size_distribution *sd;
	double operator() (double radius) const { return (*sd)(radius)*radius*radius*radius; }
};

class dust_component
{
protected:
	static const int z_lim = 10000;
	int nb_ph_en, nb_of_temp, nb_of_z, verbosity;
	int *ch_arr;
	
	double t;  // auxiliary variable to store temperature, is necessary in operator();
	// energy stored in arrays has dimension [erg];
	double *ph_en_arr, *temp_arr, *abs_perH_coeff, *abs_perH_coeff_deriv, *ext_perH_coeff, *ext_perH_coeff_deriv, 
		*abs_coeff, *abs_coeff_deriv, *int_emiss, *int_emiss_deriv, *isrf_uv_phel_rate, *isrf_uv_phel_energy, 
		*isrf_vis_phel_rate, *isrf_vis_phel_energy, *cr_phel_rate, *cr_phel_energy, *isrf_uv_phel_rate_deriv, *isrf_uv_phel_energy_deriv, 
		*isrf_vis_phel_rate_deriv, *isrf_vis_phel_energy_deriv, *cr_phel_rate_deriv, *cr_phel_energy_deriv;

	void calc_int_emiss();
	void calc_derivatives();
	void del_phel_data();

public:
	int zmin, zmax;
	// wvl_exp is the exponent in the opacity dependence on the wavelength at large wavelengths;
	// area is pi*a*a, mass of one grain in g
	// heat capacity constant must be defined (debye model of heat capacity of grains, 
	// Draine & Li, ApJ 551, p.807, 2001; Cuppen et al., MNRAS 367, p.1757, 2006)
	// Cv = 12*M_PI^4/5 N k_B/(T_D)^3 *T^3, N - number of atoms in the grain;
	double radius, area, mass, conc_perH, area_perH, volume_perH, dust_mass_perH, density, wvl_exp, 
		heat_capacity_const;
	
	std::string	name;	
	void set_name(std::string n) { name = n; }
	void calc_phel_emiss(photoel_emission *, double cr_uv_flux); // min and max charge values are calculated here;

	// functions calculates the dust emissivity, energy in cm-1, temperature in K;
	// per H nuclei and divided by 2hv/lambda^2, in [cm2]:
	double emissivity_perH(double energy, double temperature) const;
	// per one grain and divided by 2hv/lambda^2, in [cm2]:
	double emissivity(double energy, double temperature) const;
	// per one grain, multiplied by energy^3:
	double operator() (double energy) const;	

	// functions return absorption and opacity per H nuclei, in [cm2]; 
	// energy in cm-1, opacity by definition is extinction;
	double absorption_perH(double energy) const;
	double opacity_perH(double energy) const; // opacity = extinction = absorption + scattering
	
	// cross section for absorption of one grain, in [cm2]:
	double absorption(double energy) const;

	// the amount of energy, that one dust grain emits [erg s-1]:
	double get_int_emiss(double temperature) const;
	
	// returned parameters: reaction rate in s-1 (rate of electron emission by one grain), electron energy in erg:
	double get_isrf_uv_phel_rate(double charge, double & energy) const;
	double get_isrf_vis_phel_rate(double charge, double & energy) const;
	double get_cr_phel_rate(double charge, double & energy) const;

	// returns nb_of_z:
	int get_nb_of_z() { return nb_of_z; };
	// returns the value of array ch_arr:
	int get_z(int i);

	dust_component(const std::string fname, int verbosity = 1);
	// Note: the ranges of size distribution must be inside the ranges of file data:
	dust_component(const std::string fname, size_distribution *, std::string name = "", int verbosity = 1);
	~dust_component();
};

class dust_model
{
public:
	int nb_of_comp, verbosity;
	double dust_mass_perH, area_perH, conc_perH;
	
	std::string name;
	std::vector<dust_component *> components;
	
	void add_component(dust_component *dc);
	void save_data(const std::string &path) const;	
	
	double absorption_perH(double energy) const;
	double opacity_perH(double energy) const;
	
	// divided by 2hv/lambda^2, the array containing temperature values of dust components must be given, [cm2]: 
	double emissivity_perH(double energy, double *temperature) const;
	double emissivity_perH(double energy, double temperature) const; // one temperature for all components;

	// divided by 2hv/lambda^2, the arrays containing grain temperatures and concentrations must be given, [cm-1]:
	double emissivity(double energy, double *temperature, double *concentration) const;
	double emissivity(double energy, const std::vector<double> & temperature, const std::vector<double> & concentration) const;
	
	// concentration of grains of different dust components must be given, [cm-1]:
	double absorption(double energy, double *concentration) const;
	double absorption(double energy, const std::vector<double> & concentration) const;
	
	// dust component parameters, there is no check if i lies in the correct range:
	double get_grain_radius(int i) const { return components[i]->radius; }
	double get_grain_area(int i) const { return components[i]->area; }
	double get_grain_mass(int i) const { return components[i]->mass; }
	
	dust_model(int verbosity);
	virtual ~dust_model();
};

// simple model, only two components: PAH molecules and large carbon (or silicate) grains;
class two_component_dust_model : public dust_model
{
public:
	// carbon nuclei abundance in PAH molecules, per H, typical value 10^{-5};
	// dust gas mass ratio is for large grains; here, pah mass density in gas is assumed to be negligibly small;
	two_component_dust_model(const std::string &path, double c_abund_pah, double dg_ratio, double he_to_h_nb_ratio, 
		double cr_uv_flux, int verbosity = 1);
};

// Draine, ARAA 41, p. 241 (2003);
class dust_synthetic_Draine2003 : public dust_model
{
public:
	// file name with path;
	dust_synthetic_Draine2003(const std::string &fname, int verbosity = 1);
};

class dust_Draine2003 : public dust_model
{
public:
	int nb_of_carb_comp, nb_of_silic_comp;
	// nb of the dust size distribution is provided, must be less than NB_OF_DUST_MODELS_WD:
	dust_Draine2003(const std::string &path, int sd_nb, double cr_uv_flux, int verbosity = 1);
};


// The radiation of the dust-gas clump;
class dust_clump_radiation : public radiation_field
{
private:
	const dust_model *dust;

public:
	double h_col_dens; // H nuclei column density;
	double get_intensity(double) const;
	dust_clump_radiation(const dust_model *, double dust_temp, double h_col_dens, double dilution);
};

// Methods of this class return the dust heating rate due to absorption of interstellar UV, VIS, radio radiation fields;
// Check the extinction law used to calculate IS radiation extinction;
class dust_heating_ISRF
{
protected:
	int nb_of_comp, nb_ve; // nb of dust components, nb of visual extinction values;
	double *extinction, *heating_cr, **heating_ra, **heating_ir, **heating_vis, **heating_uv, **heating_fuv; // extinction in magnitudes;
	std::string name;

public:
	// calculates the table of heating rates for dust model given (dust may have several components);
	// the synthetic extinction curves are used in calculating of IS radiation attenuation (Draine, ARAA 41, p. 241, 2003);
	void calculate(const std::string &path, const dust_model *, double cr_uv_flux);

	// returns dust heating rate for each dust component, based on visiual extinction, uv field strength, cr ioniz rate;
	// heating rate units are erg s-1, heating rate of one grain, number of dust components must be given;
	void get_heating(double visual_extinct, double ir_field_strength, double uv_field_strength, double cr_ioniz_rate, 
		double *heating_rate, int nb) const;
	void save_data(const std::string &path) const;

	dust_heating_ISRF();
	~dust_heating_ISRF();
};

// Auxulary class needed for the calculation of dust heating rate;
class rf_func
{
private:
	const radiation_field *rfield;
	const dust_model *dust, *extinct_law;

public:
	int i; // dust component index
	double l0; //
	
	double operator() (double energy) const;
	void set(const radiation_field *rf, const dust_model *d, const dust_model *extl) { rfield = rf; dust = d; extinct_law = extl; }
	rf_func() : rfield(0), dust(0), extinct_law(0), l0(0.), i(0) {;}
};
