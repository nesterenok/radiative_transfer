// Changes:
// 1. The interstellar radiation field was added;
// 2. class of dust radiation was removed to the dust_model.cpp file;
// 3. check for errors 02.03.2017;
// 4. check for errors 29.08.2017;
// 5. check for errors 30.01.2017;
// 6. check for errors 14.03.2018; dust emission was added to the ISRF_Mathis1983 class;

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <fstream>
#include <iostream>
#include <cfloat>
#include "utils.h"
#include "radiation_field.h"
#include "constants.h"
#include "integration.h"

#define SOURCE_NAME "radiation_field.cpp"
using namespace std;

// some arbitrary large values for minimal and maximal energies of radiation are set;
radiation_field::radiation_field() : temperature(0.), lim_cos(0.), name(""), en_min(1.e-8), en_max(1.e+8)
{;}

double radiation_field::get_intensity(double energy, double cos_theta) const
{
	if (cos_theta < lim_cos) 
		return 0.;
	else return get_intensity(energy);
}

// emin, emax in cm-1
void radiation_field::print_parameters(double emin, double emax) const
{
	double a;
	radiation_energy_density red(this);
	radiation_photon_flux rpf(this);

	cout << "Parameters of the radiation field: " << name << endl;

	// calculation of angle- and frequency-integrated energy density, in erg cm-3;
	// G0 = 1 corresponds to 5.29e-14 erg cm-3, E = 6-13.6 eV (Draine 2011, p.123);
	a = 8.*M_PI*PLANCK_CONSTANT*SPEED_OF_LIGHT 
		*qromb<radiation_energy_density>(red, emin, emax, 1.e-6);	
	cout << "Energy density at " << emin << " - " << emax << " cm-1, in erg cm-3: " << a << endl;
	
	// calculation of angle- and frequency-integrated photon number flux, in ph cm-2 s-1;
	a = 8.*M_PI*SPEED_OF_LIGHT 
		*qromb<radiation_photon_flux>(rpf, emin, emax, 1.e-6);
	cout << "Photon number flux at " << emin << " - " << emax << " cm-1, ph cm-2 s-1: " << a << endl << endl;
}

// Radiation intensity is normalized by 2hv/lambda^2; energy is in cm-1, temperature in K;
double blackbody_radiation_field::get_intensity(double energy) const
{
	return 1./(exp(energy *CM_INVERSE_TO_KELVINS/temperature) - 1.);
}

stellar_radiation_field::stellar_radiation_field(double t, double r, double d)
: blackbody_radiation_field(t), star_radius(r), distance(d)
{
	lim_cos = sqrt(distance*distance - star_radius *star_radius) /distance;
}

//
// Interstellar radiation field models
//
ISRF_Mathis1983::ISRF_Mathis1983() 
{
	name = "ISRF(Mathis_A&A_v128_p212_1983)";
	en_max = 40816.; // 5.06 eV = 2450 A = 40816.3 cm-1
}

double ISRF_Mathis1983::get_intensity(double energy) const
{
	if (energy < en_max) {
		double a = 7.e-13/(exp(energy *CM_INVERSE_TO_KELVINS/3000.) - 1.) + 1.65e-13/(exp(energy *CM_INVERSE_TO_KELVINS/4000.) - 1.)
			+ 1.e-14/(exp(energy *CM_INVERSE_TO_KELVINS/7500.) - 1.);
		// CMB
		a += 1./(exp(energy *CM_INVERSE_TO_KELVINS/2.728) - 1.);
		// Dust emission (Hocuk et al., 2017)
		a += 3.4e-9/(exp(energy *CM_INVERSE_TO_KELVINS/250.) - 1.) + 2.e-4/(exp(energy *CM_INVERSE_TO_KELVINS/23.3) - 1.);
		return a;
	}
	return 0.;
}

ISRF_UV_Mathis1983::ISRF_UV_Mathis1983()
{
	double c = 1./(8.*M_PI*PLANCK_CONSTANT*SPEED_OF_LIGHT);
	name = "ISRF_UV(Mathis_A&A_v128_p212_1983)";
	
	a1 = 2.373e-14 *pow(10., -4*0.6678)*c; // 1 mu = 1.e-4 cm = 10000 cm-1
	a2 = 6.825e-13 *1.e+4*c;
	a3 = 1.287e-9 *pow(10., 4*4.4172)*c;

	en_min = 40816.; // 2450 A = 40816.3 cm-1
}

// Radiation energy must be in cm-1; 1 A = 1.e-8 cm 
double ISRF_UV_Mathis1983::get_intensity(double energy) const
{
	if (energy > en_min) {
		if (energy < 74627.) { // 1340 A = 74627 cm-1
			return a1*pow(energy, -3.3322); // 0.6678 - 4 = -3.3322
		}
		else if (energy < 90909.) { // 1100 A = 90909 cm-1
			return a2*pow(energy, -5.);
		}
		else if (energy < 109650.) { // 912 A = 109649.1 cm-1
			return a3*pow(energy, -8.4172);
		}
	}
	return 0.;
}

ISRF_UV_Draine1978::ISRF_UV_Draine1978()
{
	// intensity in erg cm-2 s-1 sr-1 Hz-1 normalized by 2hv/l^2;
	// 1 [eV] = EV_TO_ERGS [erg] = EV_TO_ERGS/PLANCK_CONSTANT [Hz]  
	double c = 0.5*PLANCK_CONSTANT/EV_TO_ERGS;
	name = "ISRF_UV(Draine_ApJS_v36_p595_1978)";

	a1 = 1.658e+6 *c*CM_INVERSE_TO_EV;
	a2 = -2.152e+5*c*CM_INVERSE_TO_EV*CM_INVERSE_TO_EV;
	a3 = 6.919e+3 *c*CM_INVERSE_TO_EV*CM_INVERSE_TO_EV*CM_INVERSE_TO_EV;

	en_min = 40816.; // 2450 A = 40816.3 cm-1, note that 5 eV = 40327.7 cm-1;
}

double ISRF_UV_Draine1978::get_intensity(double energy) const
{
	if (energy > en_min && energy < 109650.) {
		return a1/energy + a2 + a3*energy;
	}
	return 0.;
}

//
// CR induced UV radiation field
//

CR_induced_UV_field::CR_induced_UV_field(double norm_factor)
{
	// 57143 < E < 117650 cm-1 or 7.09 - 14.6 eV, 
	en_min = 5.7e+4;
	en_max = 1.18e+5; 

    // Ibanez-Meja et al. 2019, E = 90330 - 109690 cm-1
    // en_min = 9.0e+4;
    // en_max = 1.1e+5;

	a = norm_factor/(8.*M_PI*(en_max - en_min) *SPEED_OF_LIGHT); 
	name = "CR_UV";
}

// energy is in cm-1;
double CR_induced_UV_field::get_intensity(double energy) const
{
	if (energy > en_min && energy < en_max) {
		return a/(energy*energy);
	}
	return 0.;
}
