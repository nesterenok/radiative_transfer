
#pragma once
#include <memory>
#include <string>

// The external radiation field with pure virtual method.
// intensity is dimensionless, I [erg cm-2 s-1 Hz-1 sr-1] normalized by 2hv/lambda^2;
class radiation_field
{
public:
	std::string name;
	double temperature, lim_cos, en_min, en_max; // temperature in K;

	// radiation energy must be in cm-1, energy > 0.:
	double operator() (double energy) const { return get_intensity(energy); }
	virtual double get_intensity(double energy) const = 0;
	virtual double get_intensity(double energy, double cos_theta) const;
	virtual double get_dilution() const { return 0.5*(1. - lim_cos); }
	
	// dilution factor must be in the range [0,1], no dilution = 1:
	virtual void set_dilution(double del) { lim_cos = 1. - 2.*del; }
	
	// print angle- and frequency-integrated energy density and photon flux, energy in cm-1:
	void print_parameters(double emin, double emax) const;

	radiation_field();
};

struct radiation_energy_density
{
	const radiation_field *rfield;
	// energy in cm-1, to get the angle-integrated energy density, one must multiply by 8*pi*planck_const*speed_light;
	// note: energy density = intensity/speed of light
	double operator() (double energy) { return (*rfield)(energy) *energy*energy*energy; }
	radiation_energy_density(const radiation_field *rf) : rfield(rf) {;}
};

struct radiation_photon_flux
{
	const radiation_field *rfield;
	// energy in cm-1, to get angle-integrated photon flux, one must multiply by 8*pi*speed_light;
	double operator() (double energy) { return (*rfield)(energy) *energy*energy; }
	radiation_photon_flux(const radiation_field *rf) : rfield(rf) {;}
};

// The blackbody radiation;
class blackbody_radiation_field : public radiation_field
{
public:
	double get_intensity(double energy) const;
	blackbody_radiation_field(double temp) { temperature = temp; lim_cos = -1.; name = "blackbody";}
};

// The radiation field of the star;
class stellar_radiation_field : public blackbody_radiation_field
{
public:
	double star_radius, distance; // distance to the star centre;
	stellar_radiation_field(double temperature, double star_radius, double distance);
};

// The interstellar radiation field; 
// Draine, Interstellar medium (2011); Mathis et al., A&A 128, p. 212 (1983); Hocuk et al., A&A 604, A58 (2017); 
// dust emission and CMB are included, photon energy < 5 eV;
// Note: Av is measured at 5550 A = 18000 cm-1 = 2.2 eV;
class ISRF_Mathis1983 : public radiation_field
{
public:
	double get_intensity(double energy) const;
	ISRF_Mathis1983();
};

// Far-ultraviolet, FUV;
// the spectrum has a factor G0 = 1.14 (relative to estimate by Habing, Bull. Astr. Inst. Netherlands, 19, 421, 1968);
// G0 = 1 corresponds to energy density 5.29e-14 erg cm-3; hydrogen atom ionization energy is 13.5984 eV;
class ISRF_UV_Mathis1983 : public radiation_field
{
private:
	double a1, a2, a3;

public:
	double get_intensity(double energy) const;
	ISRF_UV_Mathis1983();
};

// The UV background in the solar neighborhood according to Draine, ApJS 36, p. 595 (1978);
// the approximation by Draine (1978) of UV field is valid for 5 eV - 13.6 eV; the spectrum has a factor G0 = 1.69; 
// see the extention of the approximation for E < 6 eV (l > 200 nm) by Dishoeck & Black, ApJ 258, p.533 (1982), Heays et al. A&A 602, id.A105 (2017)
class ISRF_UV_Draine1978 : public radiation_field
{
private:
	double a1, a2, a3;

public:
	double get_intensity(double) const;
	ISRF_UV_Draine1978();
};

// Cosmic ray induced FUV radiation field;
// the photon flux is assumed to be constant at 850 - 1750 A, or 57143 < E < 117650 cm-1, or 7.09 - 14.6 eV (Gredel et al. ApJ 347, p. 289, 1989);
// it is assumed that the emission probability for photon of any energy is equal to a constant value;
// Ibanez-Meja et al. 2019, MNRAS (accepted) used the energy range for CR induced UV radiation 11.2 - 13.6 eV (90330 - 109690 cm-1)
class CR_induced_UV_field : public radiation_field
{
private:
	double a;

public:
	double get_intensity(double) const;
	CR_induced_UV_field(double);
};
