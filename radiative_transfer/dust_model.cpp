
// Be carefull, changes:
// 1. The classes of grain size distribution were added; all other classes were rewritten; 
// 2. The fucntion absorption() was added instead of opacity(), now opacity is equal to extinction;
// 3. The new class dust_sinthetic_Draine2003 was added instead of old class where synthetic curve was used;
// 10.03.2017. Check for errors. Extinction law was added in the calculations of dust heating by IS radiation;
// 01.09.2017. Check for errors.
// 27.11.2017. Error was found: magnitude = 1.086* optical depth;
// 29.01.2017. Check for errors.

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <stdlib.h>
#include <memory.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <float.h>
#include <sstream>
#include <algorithm>

#include "constants.h"
#include "utils.h"
#include "integration.h"
#include "special_functions.h"
#include "dust_model.h"

#define MAX_TEXT_LINE_WIDTH 240
#define SOURCE_NAME "dust_model.cpp"
using namespace std;

//
// Grain size distributions
//

size_distribution::size_distribution(double r) : rad(r), rad_min(r), rad_max(r), norm(0.) {
	distribution = false;
}

size_distribution::size_distribution(double rmin, double rmax) : rad(rmin), rad_min(rmin), rad_max(rmax), norm(0.) {
	distribution = true;
}

size_distribution_mono::size_distribution_mono(double r, double conc) : size_distribution(r) {
	norm = conc;
}

size_distribution_Draine2001::size_distribution_Draine2001(int model_nb, double rmin, double rmax) 
	: size_distribution(rmin, rmax), at(0.), ac(0.), alpha(0.), beta(0.) 
{
	if (rad_min < 3.5e-8) // in cm;
		rad_min = 3.5e-8;

	// Rv = 3.1 scaling factor 0.93; Rv = 4.0 factor 1.18; Rv = 5.5 factor 1.42 (Draine, ARA&A 41, p. 241, 2003);
	// http://www.astro.princeton.edu/~draine/dust/dustmix.html
	switch (model_nb) 
	{
		case 0:	scaling_factor = 0.93; break;
		case 1:
		case 2: scaling_factor = 1.42; break;
		default: scaling_factor = 1.;
	}
}

double size_distribution_Draine2001::operator() (double radius) const
{
	double a, answ = norm*pow(radius/at, alpha)/radius;
	
	if (beta < 0)
		answ /= 1. - beta*radius/at;
	else answ *= 1. + beta*radius/at;

	if (radius > at) {
		a = (radius - at)/ac;
		answ *= exp(-a*a*a);
	}
	return answ;
}

// For the parameter values see Weingartner & Draine(2001) and Draine & Li (2007); 1 A = 1.e-8 cm;
size_distribution_carbon_Draine2001::size_distribution_carbon_Draine2001(int m_nb, double rmin, double rmax)
	: size_distribution_Draine2001(m_nb, rmin, rmax), a1(4.e-8), a2(2.e-7), sigma1(0.4), sigma2(0.55)
{
	err_func ef;
	double p, density;
	
	// parameters of lognormal terms:
	carbon_abundance = scaling_factor*WD2001_model_data[m_nb][1];
	
	density = 2.24; // in g/cm3
	p = 3.*12.*ATOMIC_MASS_UNIT*carbon_abundance/(pow(2.*M_PI, 1.5)*density);

	b1 = 0.75*p*exp(-4.5*sigma1*sigma1)/(a1*a1*a1*sigma1 *( 1. + ef.f(M_SQRT1_2*(3.*sigma1 + log(a1/3.5e-8)/sigma1)) ));
	b2 = 0.25*p*exp(-4.5*sigma2*sigma2)/(a2*a2*a2*sigma2 *(1. + ef.f(M_SQRT1_2*(3.*sigma2 + log(a2/3.5e-8)/sigma2))));

	// parameters of non-lognormal terms:
	alpha = WD2001_model_data[m_nb][2];
	beta = WD2001_model_data[m_nb][3];
	
	at = 1.e-4*WD2001_model_data[m_nb][4]; // in cm;
	ac = 1.e-4*WD2001_model_data[m_nb][5];		
	norm = scaling_factor*WD2001_model_data[m_nb][6];
}

double size_distribution_carbon_Draine2001::operator() (double radius) const 
{
	double answ;
	// non-lognormal term:
	answ = size_distribution_Draine2001::operator()(radius);
	
	// if radius (in cm) is small enough, lognormal terms are calculated:
	if (radius < 3.e-6) 
	{
		double c1, c2;
		c1 = log(radius/a1)/sigma1;
		c2 = log(radius/a2)/sigma2;
		answ += (b1*exp(-0.5*c1*c1) + b2*exp(-0.5*c2*c2))/radius;
	}
	return answ;
}

size_distribution_silicate_Draine2001::size_distribution_silicate_Draine2001(int m_nb, double rmin, double rmax)
	: size_distribution_Draine2001(m_nb, rmin, rmax)
{
	alpha = WD2001_model_data[m_nb][7];
	beta = WD2001_model_data[m_nb][8];
	
	at = 1.e-4*WD2001_model_data[m_nb][9]; // in cm	
	ac = 1.e-5;			
	norm = scaling_factor*WD2001_model_data[m_nb][10];
}

//
// Dust component class
//

// Data on synthetic extinction curves are read;
dust_component::dust_component(const std::string fname, int verb) : verbosity(verb),
	ch_arr(0), isrf_uv_phel_rate(0), isrf_uv_phel_energy(0), isrf_vis_phel_rate(0), isrf_vis_phel_energy(0), 
	cr_phel_rate(0), cr_phel_energy(0), isrf_uv_phel_rate_deriv(0), isrf_uv_phel_energy_deriv(0), isrf_vis_phel_rate_deriv(0),
	isrf_vis_phel_energy_deriv(0), cr_phel_rate_deriv(0), cr_phel_energy_deriv(0), nb_of_z(0), zmin(0), zmax(0), t(1.), radius(0.), 
	mass(0.), volume_perH(0.), density(0.), heat_capacity_const(0.)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i;
	double a, th, th2, albedo, c_ext, k_abs, mdr;
	
	ifstream input;
	input.open(fname.c_str(), ios_base::in);	
	
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}
	
	// comments:
	for (i = 0; i < 6; i++) {
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	}

	input >> name >> dust_mass_perH >> mdr >> area_perH >> conc_perH >> wvl_exp >> text_line;
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_ph_en;

	if (verbosity) {
		cout << "Dust component is initializing: " << name << endl;
	}

	ph_en_arr = new double [nb_ph_en];
	abs_coeff = new double [nb_ph_en];
	abs_perH_coeff = new double [nb_ph_en];
	ext_perH_coeff = new double [nb_ph_en];
	
	for (i = 0; i < nb_ph_en; i++)
	{
		input >> a >> albedo >> th >> c_ext >> k_abs >> th2;
		ph_en_arr[i] = 1.e+4/a;						// um to cm-1 conversion;
		abs_perH_coeff[i] = k_abs *dust_mass_perH;	// absorption in the file is in cm2/g;
		ext_perH_coeff[i] = c_ext;					// cm2 per H;
		
		abs_coeff[i] = abs_perH_coeff[i]/conc_perH;  // cm2, cross section for absorption of one grain
	}
	input.close();
	
	calc_derivatives(); // this function must be called before the next one;
	calc_int_emiss();
}

dust_component::dust_component(const std::string fname, size_distribution *s_distr, string str, int verb) : verbosity(verb),
	ch_arr(0), isrf_uv_phel_rate(0), isrf_uv_phel_energy(0), isrf_vis_phel_rate(0), isrf_vis_phel_energy(0), 
	cr_phel_rate(0), cr_phel_energy(0), isrf_uv_phel_rate_deriv(0), isrf_uv_phel_energy_deriv(0), isrf_vis_phel_rate_deriv(0), 
	isrf_vis_phel_energy_deriv(0), cr_phel_rate_deriv(0), cr_phel_energy_deriv(0), nb_of_z(0), zmin(0), zmax(0), t(1.), heat_capacity_const(0.)
{
	char text_line[MAX_TEXT_LINE_WIDTH];	
	int i, j, nb_gr;
	double r1, r2, d, q1, q2, in1, in2;
	double *q_abs, *q_sca, *q_cth, *gr_rad_arr;
	
	sd_func0 f0;
	sd_func1 f1;
	sd_func2 f2;
	sd_func3 f3;

	ifstream input;
	input.open(fname.c_str(), ios_base::in);
	
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> name >> density >> wvl_exp;
	input >> nb_ph_en >> nb_gr;

	if (str != "")
		name = str;

	if (verbosity) {
		cout << "Dust component is initializing: " << name << endl;
	}
	
	ph_en_arr = new double [nb_ph_en];
	abs_perH_coeff = new double [nb_ph_en];
	ext_perH_coeff = new double [nb_ph_en];
	abs_coeff = new double [nb_ph_en];
	
	q_abs = new double [nb_gr];
	q_sca = new double [nb_gr];
	q_cth = new double [nb_gr];
	gr_rad_arr = new double [nb_gr];

	for (i = 0; i < nb_gr; i++) {
		input >> gr_rad_arr[i];
	}

	if (s_distr->distribution)
	{
		if (verbosity && (gr_rad_arr[0] > s_distr->rad_min || gr_rad_arr[nb_gr-1] < s_distr->rad_max)) {
			cout << "  Warning: size distribution is out of range of file data: " << name << endl;
		}
		f0.sd = f1.sd = f2.sd = f3.sd = s_distr;
	
		conc_perH   = qromb<sd_func0>(f0, s_distr->rad_min, s_distr->rad_max, 1.e-6);
		radius      = qromb<sd_func1>(f1, s_distr->rad_min, s_distr->rad_max, 1.e-6);
		area_perH   = qromb<sd_func2>(f2, s_distr->rad_min, s_distr->rad_max, 1.e-6) *M_PI;       // pi*a*a per H nuclei of gas;
		volume_perH = qromb<sd_func3>(f3, s_distr->rad_min, s_distr->rad_max, 1.e-6) *4./3.*M_PI; // volume per H nuclei;
	
		dust_mass_perH = density*volume_perH;
	
		// radius, area and mass of one grain:
		radius = radius/conc_perH;
		area = area_perH/conc_perH;
		mass = dust_mass_perH/conc_perH;

		// concentration is recalculated in order to conserve total dust volume:
		// conc_perH = 3.*volume_perH/(4.*M_PI *radius*radius*radius);
	
		i = 0;
		while (!input.eof() && i < nb_ph_en)
		{
			q1 = q2 = 0.;
			input >> ph_en_arr[i];
			input >> q_abs[0] >> q_sca[0] >> q_cth[0];

			for (j = 1; j < nb_gr; j++) 
			{
				input >> q_abs[j] >> q_sca[j] >> q_cth[j];
				if (gr_rad_arr[j-1] <= s_distr->rad_max && gr_rad_arr[j] > s_distr->rad_min) 
				{
					if (gr_rad_arr[j-1] < s_distr->rad_min)
						r1 = s_distr->rad_min;
					else r1 = gr_rad_arr[j-1];

					if (gr_rad_arr[j] > s_distr->rad_max)
						r2 = s_distr->rad_max;
					else r2 = gr_rad_arr[j];

					in1 = qromb<sd_func2>(f2, r1, r2, 1.e-4);
					in2 = qromb<sd_func3>(f3, r1, r2, 1.e-4);
			
					d = (q_abs[j] - q_abs[j-1])/(gr_rad_arr[j] - gr_rad_arr[j-1]); // derivative;

					q1 += (q_abs[j-1] - d*gr_rad_arr[j-1]) *in1;
					q1 += d *in2;

					d = (q_sca[j] - q_sca[j-1])/(gr_rad_arr[j] - gr_rad_arr[j-1]); // derivative;

					q2 += (q_sca[j-1] - d*gr_rad_arr[j-1]) *in1;
					q2 += d *in2;

					// this extrapolation is inaccurate:
					if (j == nb_gr-1 && gr_rad_arr[j] <= s_distr->rad_max) 
					{
						in1 = qromb<sd_func2>(f2, r2, s_distr->rad_max, 1.e-4);
						q1 += q_abs[j] *in1;
						q2 += q_sca[j] *in1;
					}
				}
			
				// this extrapolation is inaccurate:
				if (j == nb_gr-1 && s_distr->rad_min >= gr_rad_arr[nb_gr-1])
				{
					r1 = s_distr->rad_min;
					r2 = s_distr->rad_max;
					in1 = qromb<sd_func2>(f2, r1, r2, 1.e-4);
				
					q1 = q_abs[j] *in1;
					q2 = q_sca[j] *in1;
				}
			}
			abs_perH_coeff[i] = M_PI *q1;
			ext_perH_coeff[i] = M_PI *(q1 + q2); // extinction = absorption + scattering
		
			abs_coeff[i] = abs_perH_coeff[i]/conc_perH;
			i++;
		}
	}
	// the case of mono-size distribution:
	else
	{	
		if (verbosity && (gr_rad_arr[0] > s_distr->rad || gr_rad_arr[nb_gr-1] < s_distr->rad)) {
			cout << "  Warning: size distribution is out of range of file data: " << name << endl;
		}
		
		conc_perH = s_distr->norm;
		radius = s_distr->rad;
		
		area = M_PI *radius*radius;
		mass = 4./3.*M_PI*radius*radius*radius*density;

		area_perH = area *conc_perH;
		volume_perH = 4./3.*M_PI *radius*radius*radius *conc_perH;
		dust_mass_perH = density*volume_perH;

		i = 0;
		while (!input.eof() && i < nb_ph_en)
		{
			input >> ph_en_arr[i];
			input >> q_abs[0] >> q_sca[0] >> q_cth[0];

			for (j = 1; j < nb_gr; j++) 
			{
				input >> q_abs[j] >> q_sca[j] >> q_cth[j];
				if (gr_rad_arr[j-1] <= radius && gr_rad_arr[j] > radius) 
				{
					q1 = q_abs[j-1] + (radius - gr_rad_arr[j-1]) *(q_abs[j] - q_abs[j-1])/(gr_rad_arr[j] - gr_rad_arr[j-1]);
					q2 = q_sca[j-1] + (radius - gr_rad_arr[j-1]) *(q_sca[j] - q_sca[j-1])/(gr_rad_arr[j] - gr_rad_arr[j-1]);
				}
				else if (j == nb_gr-1 && radius >= gr_rad_arr[j])
				{
					q1 = q_abs[j];
					q2 = q_sca[j];
				}
				else if (j == 1 && radius < gr_rad_arr[j-1])
				{
					q1 = q_abs[0];
					q2 = q_sca[0];
				}
			}

			abs_perH_coeff[i] = area_perH *q1;
			ext_perH_coeff[i] = area_perH *(q2 + q1);
		
			abs_coeff[i] = area *q1; // per one grain, [cm2]
			i++;
		}
	}
	input.close();
	
	calc_derivatives(); // this function must be called before the next one;
	calc_int_emiss();
	
	delete [] q_abs;
	delete [] q_sca;
	delete [] q_cth;
	delete [] gr_rad_arr;
}

dust_component::~dust_component()
{
	delete [] ph_en_arr;
	delete [] temp_arr;
	delete [] ch_arr;
	
	delete [] abs_coeff; 
	delete [] abs_perH_coeff;
	delete [] ext_perH_coeff;
	delete [] int_emiss;

	delete [] abs_coeff_deriv;
	delete [] abs_perH_coeff_deriv;
	delete [] int_emiss_deriv;
	delete [] ext_perH_coeff_deriv;

	del_phel_data();
}

void dust_component::del_phel_data()
{
	delete [] isrf_uv_phel_rate;
	delete [] isrf_uv_phel_energy;
	delete [] isrf_vis_phel_rate;
	delete [] isrf_vis_phel_energy;
	delete [] cr_phel_rate;
	delete [] cr_phel_energy;

	delete [] isrf_uv_phel_rate_deriv; 
	delete [] isrf_uv_phel_energy_deriv; 
	delete [] isrf_vis_phel_rate_deriv; 
	delete [] isrf_vis_phel_energy_deriv; 
	delete [] cr_phel_rate_deriv; 
	delete [] cr_phel_energy_deriv;
}

double dust_component::emissivity_perH(double energy, double temp) const {
	return absorption_perH(energy)/(exp(energy *CM_INVERSE_TO_KELVINS/temp) - 1.);
}

// t - is internal variable, must be initialized before function call; 
double dust_component::operator() (double energy) const {
	return absorption(energy)/(exp(energy *CM_INVERSE_TO_KELVINS/t) - 1.) *energy*energy*energy;
}

double dust_component::emissivity(double energy, double temp) const {
	return absorption(energy)/(exp(energy *CM_INVERSE_TO_KELVINS/temp) - 1.);
}

double dust_component::absorption_perH(double energy) const
{
	if (energy < ph_en_arr[0])
		return abs_perH_coeff[0] *pow(ph_en_arr[0]/energy, wvl_exp);
	
	if (energy > ph_en_arr[nb_ph_en-1])
		return abs_perH_coeff[nb_ph_en-1];

	int i, l = 0, r = nb_ph_en-1;
	while (r-l > 1)
	{
		i = l + ((r - l) >> 1);
		if (ph_en_arr[i] < energy) 
			l = i;
		else r = i;
	}
	return abs_perH_coeff[l] + abs_perH_coeff_deriv[l] *(energy - ph_en_arr[l]);
}

double dust_component::opacity_perH(double energy) const
{
	if (energy < ph_en_arr[0])
		return ext_perH_coeff[0] *pow(ph_en_arr[0]/energy, wvl_exp);
	
	if (energy > ph_en_arr[nb_ph_en-1])
		return ext_perH_coeff[nb_ph_en-1];

	int i, l = 0, r = nb_ph_en-1;
	while (r-l > 1)
	{
		i = l + ((r - l) >> 1);
		if (ph_en_arr[i] < energy) 
			l = i;
		else r = i;
	}
	return ext_perH_coeff[l] + ext_perH_coeff_deriv[l] *(energy - ph_en_arr[l]);
}

double dust_component::absorption(double energy) const
{
	if (energy < ph_en_arr[0])
		return abs_coeff[0] *pow(ph_en_arr[0]/energy, wvl_exp);
	
	if (energy > ph_en_arr[nb_ph_en-1])
		return abs_coeff[nb_ph_en-1];

	int i, l = 0, r = nb_ph_en-1;
	while (r-l > 1)
	{
		i = l + ((r - l) >> 1);
		if (ph_en_arr[i] < energy) 
			l = i;
		else r = i;
	}
	return abs_coeff[l] + abs_coeff_deriv[l] *(energy - ph_en_arr[l]);
}

double dust_component::get_int_emiss(double temp) const
{
	if (temp < temp_arr[0])
		return int_emiss[0];

	if (temp > temp_arr[nb_of_temp-1])
		return int_emiss[nb_of_temp-1];

	int i, l = 0, r = nb_of_temp-1;
	while (r-l > 1)
	{
		i = l + ((r - l) >> 1);
		if (temp_arr[i] < temp) 
			l = i;
		else r = i;
	}
	return int_emiss[l] + int_emiss_deriv[l] *(temp - temp_arr[l]);
}

// All necessary data must be initialized at this point:
void dust_component::calc_derivatives()
{
	abs_coeff_deriv = new double [nb_ph_en-1];
	abs_perH_coeff_deriv = new double [nb_ph_en-1];
	ext_perH_coeff_deriv = new double [nb_ph_en-1];
	
	for (int i = 0; i < nb_ph_en-1; i++) 
	{
		abs_coeff_deriv[i] = (abs_coeff[i+1] - abs_coeff[i])/(ph_en_arr[i+1] - ph_en_arr[i]);
		abs_perH_coeff_deriv[i] = (abs_perH_coeff[i+1] - abs_perH_coeff[i])/(ph_en_arr[i+1] - ph_en_arr[i]);
		ext_perH_coeff_deriv[i] = (ext_perH_coeff[i+1] - ext_perH_coeff[i])/(ph_en_arr[i+1] - ph_en_arr[i]);
	}
}

void dust_component::calc_int_emiss()
{
	int i, nb_t;
	double en_min, t_min, t_max, t_sc;

	nb_t = 100;
	t_sc = pow(10., 1./nb_t);
	
	t_min = 1.;
	t_max = 3000.;
	nb_of_temp = (int) (nb_t *log10(t_max/t_min));
	
	temp_arr = new double [nb_of_temp];
	int_emiss = new double [nb_of_temp];

	if (verbosity) {
		cout << "  calculation of integrated emissivity for dust component " << name << " ..." << endl;
	}
	
	en_min = 0.01;	// in cm^{-1}
	for (i = 0, t = t_min; i < nb_of_temp; i++, t *= t_sc)
	{
		temp_arr[i] = t;
		// the upper limit was chosen to be flexible:
		int_emiss[i] = qromb<dust_component>(*this, en_min, 40.*t, 1.e-6) 
			*8*M_PI*PLANCK_CONSTANT *SPEED_OF_LIGHT *SPEED_OF_LIGHT;
	}
	
	int_emiss_deriv = new double [nb_of_temp-1];
	for (i = 0; i < nb_of_temp-1; i++) {
		int_emiss_deriv[i] = (int_emiss[i+1] - int_emiss[i])/(temp_arr[i+1] - temp_arr[i]);
	}
}

void dust_component::calc_phel_emiss(photoel_emission *ph_emiss, double cr_uv_flux)
{
	int i, j;
	double rate, energy;
	vector<int> chs1, chs2;

	// interstellar radiation field, E > 40816.3 cm-1 (2450 A, 5.06 eV);
	ISRF_UV_Draine1978 *rf1 
		= new ISRF_UV_Draine1978();

	// interstellar radiation field, E < 40816.3 cm-1;
	ISRF_Mathis1983 *rf2
		= new ISRF_Mathis1983();

	CR_induced_UV_field *rf3 
		= new CR_induced_UV_field(cr_uv_flux);
	
	if (verbosity) {
		cout << "  calculation of photoelectron emission rate for dust component " << name << " ..." << endl;
	}
	
	i = 0;
	chs1.push_back(0);
	chs2.push_back(0);
	ph_emiss->get_limit_charges(radius, zmin, zmax);

	do {
		i++;
		if (i <= 30)
			j = i;
		else if (i <= 100)
			j = 30 + 2*(i - 30);
		else if (i <= 300)
			j = 170 + 5*(i - 100);
		else if (i <= 500)
			j = 1170 + 10*(i - 300);
		else j = 3170 + 50*(i - 500);

		if (j <= zmax)
			chs1.push_back(j);

		if (-j >= zmin)
			chs2.push_back(-j);
	}
	while ((j < zmax || -j > zmin) && j < z_lim);

	reverse(chs2.begin(), chs2.end());
	nb_of_z = (int) chs1.size() + (int) chs2.size() - 1;
	
	delete [] ch_arr;
	ch_arr = new int [nb_of_z];
	
	for (i = 0; i < (int) chs2.size(); i++) {
		ch_arr[i] = chs2[i];
	}
	for (i = 1; i < (int) chs1.size(); i++) {
		ch_arr[(int) chs2.size() + i - 1] = chs1[i];
	}

	// arrays are deleted first,
	del_phel_data();

	isrf_uv_phel_rate = new double [nb_of_z];
	isrf_uv_phel_energy = new double [nb_of_z];
	isrf_vis_phel_rate = new double [nb_of_z];
	isrf_vis_phel_energy = new double [nb_of_z];
	cr_phel_rate = new double [nb_of_z];
	cr_phel_energy = new double [nb_of_z];

	ph_emiss->set_rfield(rf1);
	for (i = 0; i < nb_of_z; i++)
	{
		ph_emiss->get_phel_emiss(ch_arr[i], radius, rate, energy); 
		isrf_uv_phel_rate[i] = rate; // [s-1], rate of emission per one grain;
		isrf_uv_phel_energy[i] = energy*CM_INVERSE_TO_ERG;
	}

	ph_emiss->set_rfield(rf2);
	for (i = 0; i < nb_of_z; i++)
	{
		ph_emiss->get_phel_emiss(ch_arr[i], radius, rate, energy); 
		isrf_vis_phel_rate[i] = rate;
		isrf_vis_phel_energy[i] = energy*CM_INVERSE_TO_ERG;
	}

	ph_emiss->set_rfield(rf3);
	for (i = 0; i < nb_of_z; i++)
	{	
		ph_emiss->get_phel_emiss(ch_arr[i], radius, rate, energy); 
		cr_phel_rate[i] = rate;
		cr_phel_energy[i] = energy*CM_INVERSE_TO_ERG;
	}

	// linear interpolation is employed: 
	isrf_uv_phel_rate_deriv = new double [nb_of_z-1];
	isrf_uv_phel_energy_deriv = new double [nb_of_z-1];
	isrf_vis_phel_rate_deriv = new double [nb_of_z-1];
	isrf_vis_phel_energy_deriv = new double [nb_of_z-1];
	cr_phel_rate_deriv = new double [nb_of_z-1];
	cr_phel_energy_deriv = new double [nb_of_z-1];

	for (i = 0; i < nb_of_z-1; i++)
	{
		isrf_uv_phel_rate_deriv[i] = (isrf_uv_phel_rate[i+1] - isrf_uv_phel_rate[i])/(ch_arr[i+1] - ch_arr[i]);
		isrf_uv_phel_energy_deriv[i] = (isrf_uv_phel_energy[i+1] - isrf_uv_phel_energy[i])/(ch_arr[i+1] - ch_arr[i]);

		isrf_vis_phel_rate_deriv[i] = (isrf_vis_phel_rate[i+1] - isrf_vis_phel_rate[i])/(ch_arr[i+1] - ch_arr[i]);
		isrf_vis_phel_energy_deriv[i] = (isrf_vis_phel_energy[i+1] - isrf_vis_phel_energy[i])/(ch_arr[i+1] - ch_arr[i]);

		cr_phel_rate_deriv[i] = (cr_phel_rate[i+1] - cr_phel_rate[i])/(ch_arr[i+1] - ch_arr[i]);
		cr_phel_energy_deriv[i] = (cr_phel_energy[i+1] - cr_phel_energy[i])/(ch_arr[i+1] - ch_arr[i]);
	}
	delete rf1;
	delete rf2;
	delete rf3;
}

double dust_component::get_isrf_uv_phel_rate(double charge, double & energy) const
{
	if (charge < ch_arr[0])
	{
		energy = isrf_uv_phel_energy[0];
		return isrf_uv_phel_rate[0];
	}
	else if (charge > ch_arr[nb_of_z-1])
	{
		energy = isrf_uv_phel_energy[nb_of_z-1];
		return isrf_uv_phel_rate[nb_of_z-1];
	}

	int j, l = 0, r = nb_of_z-1; 
	while (r-l > 1)
	{
		j = l + ((r-l) >> 1);
		if (ch_arr[j] < charge) 
			l = j;
		else r = j;
	}
	
	energy = isrf_uv_phel_energy[l] + isrf_uv_phel_energy_deriv[l] *(charge - ch_arr[l]);
	return isrf_uv_phel_rate[l] +  isrf_uv_phel_rate_deriv[l] *(charge - ch_arr[l]);
}

double dust_component::get_isrf_vis_phel_rate(double charge, double & energy) const
{
	if (charge < ch_arr[0])
	{
		energy = isrf_vis_phel_energy[0];
		return isrf_vis_phel_rate[0];
	}
	else if (charge > ch_arr[nb_of_z-1])
	{
		energy = isrf_vis_phel_energy[nb_of_z-1];
		return isrf_vis_phel_rate[nb_of_z-1];
	}

	int j, l = 0, r = nb_of_z-1; 
	while (r-l > 1)
	{
		j = l + ((r-l) >> 1);
		if (ch_arr[j] < charge) 
			l = j;
		else r = j;
	}
	
	energy = isrf_vis_phel_energy[l] + isrf_vis_phel_energy_deriv[l] *(charge - ch_arr[l]);
	return isrf_vis_phel_rate[l] +  isrf_vis_phel_rate_deriv[l] *(charge - ch_arr[l]);
}

double dust_component::get_cr_phel_rate(double charge, double & energy) const
{
	if (charge < ch_arr[0])
	{
		energy = cr_phel_energy[0];
		return cr_phel_rate[0];
	}
	else if (charge > ch_arr[nb_of_z-1])
	{
		energy = cr_phel_energy[nb_of_z-1];
		return cr_phel_rate[nb_of_z-1];
	}
	
	int j, l = 0, r = nb_of_z-1; 
	while (r-l > 1)
	{
		j = l + ((r-l) >> 1);
		if (ch_arr[j] < charge) 
			l = j;
		else r = j;
	}
	
	energy = cr_phel_energy[l] +  cr_phel_energy_deriv[l] *(charge - ch_arr[l]);
	return cr_phel_rate[l] +  cr_phel_rate_deriv[l] *(charge - ch_arr[l]);
}

int dust_component::get_z(int i)
{ 
	if (i >= 0 && i < nb_of_z) return ch_arr[i]; 
	else return 0; 
};

//
// Dust model class
//

dust_model::dust_model(int verb) : nb_of_comp(0), dust_mass_perH(0.), area_perH(0.), conc_perH(0.), verbosity(verb)
{;}

void dust_model::add_component(dust_component *dc) 
{ 
	components.push_back(dc);
	dust_mass_perH += dc->dust_mass_perH;
	area_perH += dc->area_perH;
	conc_perH += dc->conc_perH;
	nb_of_comp++;

	if (verbosity) {
		cout << "  dust component is added to the model: " << dc->name << endl;
	}
}

dust_model::~dust_model() 
{
	for (int i = 0; i < (int) components.size(); i++) {
		delete components[i];
	}
}

double dust_model::absorption_perH(double energy) const
{
	double a(0.);
	for (int i = 0; i < (int) components.size(); i++) {
		a += components[i]->absorption_perH(energy);
	}
	return a;
}

double dust_model::opacity_perH(double energy) const
{
	double a(0.);
	for (int i = 0; i < (int) components.size(); i++) {
		a += components[i]->opacity_perH(energy);
	}
	return a;
}

double dust_model::emissivity_perH(double energy, double *temperature) const 
{
	double a(0.);
	for (int i = 0; i < (int) components.size(); i++) {
		a += components[i]->emissivity_perH(energy, temperature[i]);
	}
	return a;
}

double dust_model::emissivity_perH(double energy, double temperature) const 
{
	double a(0.);
	for (int i = 0; i < (int) components.size(); i++) {
		a += components[i]->emissivity_perH(energy, temperature);
	}
	return a;
}

double dust_model::absorption(double energy, double *concentration) const
{
	double a(0.);
	for (int i = 0; i < (int) components.size(); i++) {
		a += components[i]->absorption(energy) *concentration[i];
	}
	return a;
}

double dust_model::absorption(double energy, const vector<double> & concentration) const
{
	double a(0.);
	for (int i = 0; i < (int)components.size(); i++) {
		a += components[i]->absorption(energy) * concentration[i];
	}
	return a;
}

double dust_model::emissivity(double energy, double *temperature, double *concentration) const
{
	double a(0.);
	for (int i = 0; i < (int) components.size(); i++) {
		a += components[i]->emissivity(energy, temperature[i]) *concentration[i];
	}
	return a;
}

double dust_model::emissivity(double energy, const vector<double> & temperature, const vector<double> & concentration) const
{
	double a(0.);
	for (int i = 0; i < (int)components.size(); i++) {
		a += components[i]->emissivity(energy, temperature[i]) * concentration[i];
	}
	return a;
}

void dust_model::save_data(const string &path) const
{
	int i, j, k;
	double a, b, sc, tot;
	string fname;
	ofstream output;
	
	fname = path + "dust_";
	fname += name;
	fname += ".txt";
	
	output.open(fname.c_str(), ios_base::out);
	output << scientific;
	output.precision(4);
	
	output << "Total area (cm2 per H): " << area_perH << endl
		<< "Dust mass (g per H): " << dust_mass_perH << endl
		<< "Grain concentration (particles per H): " << conc_perH << endl
		<< "Nb of components: " << nb_of_comp << endl;

	output << left << setw(5) << "nb" << setw(8) << "name" << setw(15) << "radius (cm)" << setw(15) << "grain mass(g)"
		 << setw(15) << "conc(per H)" << setw(15) << "area(per H)" << endl;

	for (i = 0; i < (int) components.size(); i++) {
		output << left << setw(5) << i << setw(8) << components[i]->name << setw(15) << components[i]->radius << setw(15) << components[i]->mass
			<< setw(15) << components[i]->conc_perH << setw(15) << components[i]->area_perH << endl;
	}
	output << endl;

	output.precision(3);
	output << "Energy (cm-1) - extinction normalized by H nuclei concentration [cm2], last column - total extinction:" << endl;
	
	sc = pow(10, 0.1);
	for (a = 1.; a < 2.1e+5; a *= sc) 
	{
		tot = 0.;
		output << left << setw(12) << a;
		
		for (i = 0; i < (int) components.size(); i++) 
		{
			b = components[i]->opacity_perH(a);
			tot += b;
			output << left << setw(12) << b;
		}
		output << left << setw(12) << tot << endl;
	}
	
	sc = pow(10, 0.02);
	output << endl << "Temperature (K) - integrated emissivity (energy emitted by one grain, erg/s):" << endl;
	for (a = 1.; a < 1.1e+3; a *= sc) 
	{
		output << left << setw(12) << a;	
		for (i = 0; i < (int) components.size(); i++) {
			output << left << setw(12) << components[i]->get_int_emiss(a);
		}
		output << endl;
	}

	output << endl << "Photoelectric rates on grains, rate of electron emission by one grain (s-1), electron energy in erg" << endl
		<< "Data: charge z; for IS UV field: rate, energy; for IS VIS field; for CR induced field;" << endl;
	for (i = 0; i < (int) components.size(); i++) 
	{
		output << components[i]->name << endl;
		for (j = 0; j < components[i]->get_nb_of_z(); j++) 
		{
			k = components[i]->get_z(j);
			if (abs(k) <= 100) 
			{
				a = components[i]->get_isrf_uv_phel_rate(k, b);
				output << left << setw(8) << k << setw(12) << a << setw(15) << b;

				a = components[i]->get_isrf_vis_phel_rate(k, b);
				output << left << setw(12) << a << setw(15) << b;

				a = components[i]->get_cr_phel_rate(k, b);
				output << left << setw(12) << a << setw(12) << b << endl;
			}
		}
	}
	output.close();
}

//
// Specific dust models
//

two_component_dust_model::two_component_dust_model(const string &path, double c_abund_pah, double dg_ratio, 
	double he_to_h_nb_ratio, double cr_uv_flux, int verb) 
	: dust_model(verb)
{
	int v;
	double pah_conc, grain_conc, atom_mass;
	string fname;

	name = "2_comp_model";
	photoel_emission_carbon_dust *phem_carb(0);
	photoel_emission_silicate_dust *phem_sil(0); 

	dust_component *dust_comp;
	size_distribution_mono *s_distr(0);

	// PAH:
	if (c_abund_pah > DBL_EPSILON)
	{
		// H/C = 0.5 in PAH molecules;
		// concentration per H, graphite density is used 2.24 g/cm3;
		pah_conc = 3.*c_abund_pah*12.5*ATOMIC_MASS_UNIT
			/(4.*2.24*M_PI *RADIUS_OF_PAH_MOLECULES*RADIUS_OF_PAH_MOLECULES*RADIUS_OF_PAH_MOLECULES); 
		
		fname = path + "dust/dust_PAHneut_Draine2001.txt";
		s_distr = new size_distribution_mono(RADIUS_OF_PAH_MOLECULES, pah_conc);
		dust_comp = new dust_component(fname, s_distr, "PAH", verbosity);
		
		phem_carb =	new photoel_emission_carbon_dust(path, v = 0);	
		dust_comp->calc_phel_emiss(phem_carb, cr_uv_flux);

		atom_mass = 8.2*ATOMIC_MASS_UNIT;
		dust_comp->heat_capacity_const = 2.4*M_PI*M_PI*M_PI*M_PI *BOLTZMANN_CONSTANT *dust_comp->mass
			/(atom_mass*GRAPHITE_DEBYE_TEMP*GRAPHITE_DEBYE_TEMP*GRAPHITE_DEBYE_TEMP);

		add_component(dust_comp);
		delete s_distr;
	}
	// large grains, silicate:
	if (MATERIAL_OF_LARGE_DUST_GRAINS) 
	{
		grain_conc = 3.*(1. + 4.*he_to_h_nb_ratio)*ATOMIC_MASS_UNIT*dg_ratio
			/(4.*M_PI *RADIUS_OF_LARGE_GRAINS *RADIUS_OF_LARGE_GRAINS *RADIUS_OF_LARGE_GRAINS *SILICATE_MATERIAL_DENSITY);
		s_distr = new size_distribution_mono(RADIUS_OF_LARGE_GRAINS, grain_conc);
		
		fname = path + "dust/dust_silicate_Draine2001.txt";
		//fname = path + "dust/dust_silicate_David1995.txt";
		dust_comp = new dust_component(fname, s_distr, "Sil", verbosity);
		
		phem_sil = new photoel_emission_silicate_dust(path, v = 0);
		dust_comp->calc_phel_emiss(phem_sil, cr_uv_flux);

		// olivine MgFeSiO4, average atomic number 24.6 a.m.u.
		atom_mass = 24.6*ATOMIC_MASS_UNIT;
		dust_comp->heat_capacity_const = 2.4*M_PI*M_PI*M_PI*M_PI *BOLTZMANN_CONSTANT *dust_comp->mass
			/(atom_mass*SILICATE_DEBYE_TEMP*SILICATE_DEBYE_TEMP*SILICATE_DEBYE_TEMP);
	}
	// graphite:
	else 
	{
		grain_conc = 3.*(1. + 4.*he_to_h_nb_ratio)*ATOMIC_MASS_UNIT*dg_ratio
			/(4.*M_PI *RADIUS_OF_LARGE_GRAINS *RADIUS_OF_LARGE_GRAINS *RADIUS_OF_LARGE_GRAINS *GRAPHITE_MATERIAL_DENSITY);
		s_distr = new size_distribution_mono(RADIUS_OF_LARGE_GRAINS, grain_conc);
		
		fname = path + "dust/dust_graphite_Draine1993.txt";
		dust_comp = new dust_component(fname, s_distr, "Gra", verbosity);
		
		if (phem_carb == 0)
			phem_carb =	new photoel_emission_carbon_dust(path, v = 0);
		dust_comp->calc_phel_emiss(phem_carb, cr_uv_flux);

		// pure graphite without H,
		atom_mass = 12.*ATOMIC_MASS_UNIT;
		dust_comp->heat_capacity_const = 2.4*M_PI*M_PI*M_PI*M_PI *BOLTZMANN_CONSTANT *dust_comp->mass
			/(atom_mass*GRAPHITE_DEBYE_TEMP*GRAPHITE_DEBYE_TEMP*GRAPHITE_DEBYE_TEMP);
	}
	add_component(dust_comp);

	delete s_distr;
	delete phem_carb;
	delete phem_sil;
}

dust_synthetic_Draine2003::dust_synthetic_Draine2003(const string &fname, int verb) : dust_model(verb)
{
	if (verbosity) {
		cout << "Dust model is initializing: " << fname << endl;
	}
	add_component(new dust_component(fname, verbosity));
	name = components[0]->name;
}

// Large grains of silicate and carbon provide a main contribution to the opacity at large wavelength (E < 10 cm-1);
dust_Draine2003::dust_Draine2003(const string &path, int sd_nb, double cr_uv_flux, int verb) : dust_model(verb)
{
	int i, nb, nb_of_bins, v;
	double sc, grain_minsize, grain_maxsize, atom_mass;
	double *radii;

	string fname;
	stringstream dc_name;

	if (sd_nb >= NB_OF_DUST_MODELS_WD)
		sd_nb = NB_OF_DUST_MODELS_WD-1;

	name = WD_model_names[sd_nb];
	if (verbosity) {
		cout << "Dust model is initializing: " << name << endl;
	}

	photoel_emission_carbon_dust *phem_carb =
		new photoel_emission_carbon_dust(path, v = 0);

	photoel_emission_silicate_dust *phem_sil =
		new photoel_emission_silicate_dust(path, v = 0);

	dust_component *dust_comp;

	nb = 8;
	sc = pow(10., 1./nb);
	
	// Carbonaceous grains;
	size_distribution_carbon_Draine2001 *carbon_distr;
	
	grain_minsize = WD2001_model_data[sd_nb][11]; // in cm;
	grain_maxsize = WD2001_model_data[sd_nb][12];

	nb_of_bins = (int) (nb*log10(grain_maxsize/grain_minsize)) + 1;
	radii = new double [nb_of_bins+1];

	radii[0] = grain_minsize;
	for (i = 0; i < nb_of_bins-1; i++) {
		radii[i+1] = radii[i]*sc;
	}
	radii[nb_of_bins] = grain_maxsize;

	for (i = 0; i < nb_of_bins; i++) 
	{
		dc_name.str("");
		dc_name << "C" << i;

		if (radii[i+1] < 1.e-6) // in cm
			fname = path + "dust/dust_PAHneut_Draine2001.txt";
		else fname = path + "dust/dust_graphite_Draine1993.txt";
		
		carbon_distr = new size_distribution_carbon_Draine2001(sd_nb, radii[i], radii[i+1]); 
		dust_comp = new dust_component(fname, carbon_distr, dc_name.str(), verbosity);
		
		atom_mass = 12.*ATOMIC_MASS_UNIT;
		dust_comp->heat_capacity_const = 2.4*M_PI*M_PI*M_PI*M_PI *BOLTZMANN_CONSTANT *dust_comp->mass
			/(atom_mass*GRAPHITE_DEBYE_TEMP*GRAPHITE_DEBYE_TEMP*GRAPHITE_DEBYE_TEMP);

		dust_comp->calc_phel_emiss(phem_carb, cr_uv_flux);	
		add_component(dust_comp);

		delete carbon_distr;
	}
	nb_of_carb_comp = nb_of_bins;
	delete [] radii;

	// initialization of dust components for silicate grains;
	size_distribution_silicate_Draine2001 *silicate_distr;
	
	grain_minsize = WD2001_model_data[sd_nb][13];
	grain_maxsize = WD2001_model_data[sd_nb][14];
	
	nb_of_bins = (int) (nb*log10(grain_maxsize/grain_minsize)) + 1;
	radii = new double [nb_of_bins+1];

	radii[0] = grain_minsize;
	for (i = 0; i < nb_of_bins-1; i++) {
		radii[i+1] = radii[i]*sc;
	}
	radii[nb_of_bins] = grain_maxsize;

	fname = path + "dust/dust_silicate_Draine2001.txt";
	for (i = 0; i < nb_of_bins; i++) 
	{
		dc_name.str("");
		dc_name << "S" << i;

		silicate_distr = new size_distribution_silicate_Draine2001(sd_nb, radii[i], radii[i+1]);
		dust_comp = new dust_component(fname, silicate_distr, dc_name.str(), verbosity);
	
		atom_mass = 24.6*ATOMIC_MASS_UNIT;
		dust_comp->heat_capacity_const = 2.4*M_PI*M_PI*M_PI*M_PI *BOLTZMANN_CONSTANT *dust_comp->mass
			/(atom_mass*SILICATE_DEBYE_TEMP*SILICATE_DEBYE_TEMP*SILICATE_DEBYE_TEMP);

		dust_comp->calc_phel_emiss(phem_sil, cr_uv_flux);
		add_component(dust_comp);
		
		delete silicate_distr;
	}
	nb_of_silic_comp = nb_of_bins;
	
	delete [] radii;
	delete phem_carb;
	delete phem_sil;
}

//
// Radiation of dust layer
//

dust_clump_radiation::dust_clump_radiation(const dust_model *du, double t, double h_cd, double dilution)
{
	name = "dust_radiation";
	dust = du;
	temperature = t;
	h_col_dens = h_cd;
	lim_cos = 1. - 2.*dilution;
}

// The light scattering is not taken into account;
double dust_clump_radiation::get_intensity(double energy) const {
	return ( 1. - exp(-h_col_dens *dust->opacity_perH(energy)) )/(exp(energy *CM_INVERSE_TO_KELVINS/temperature) - 1.);
}

//
// Class calculating dust heating 
//

dust_heating_ISRF::dust_heating_ISRF() : nb_of_comp(0), heating_ra(0), heating_ir(0), heating_vis(0), heating_uv(0), 
	heating_fuv(0), heating_cr(0), name("")
{
	int i, nb_d;
	double sc;
	
	nb_d = 30;
	nb_ve = (int) (3.*nb_d); // nb of visual extinction values;
	sc = pow(10., 1./nb_d);

	extinction = new double [nb_ve];

	extinction[0] = 0.;
	extinction[1] = 0.1;
	for (i = 2; i < nb_ve; i++) {
		extinction[i] = extinction[i-1]*sc;
	}
}

// At small visual extinction Av < 0.6 and G0 = 1.69, dust heating (PAH) is mainly due to UV photons (> 50%), 
// at larger Av, dust heating due to diluted starlight dominates, dust heating due to CMB is negligible, < 0.5% at Av = 20;
void dust_heating_ISRF::calculate(const std::string &path, const dust_model *dust, double cr_uv_flux)
{
	int i, j, verb;
	double a;

	string fname;
	rf_func f;
	
	name = dust->name;
	nb_of_comp = dust->nb_of_comp;
	
	if (heating_ra != 0) {
		free_2d_array(heating_ra);
		free_2d_array(heating_ir);
		free_2d_array(heating_vis);
		free_2d_array(heating_uv);
		free_2d_array(heating_fuv);
		delete [] heating_cr;
	}
	
	heating_ra = alloc_2d_array<double>(nb_of_comp, nb_ve);
	heating_ir = alloc_2d_array<double>(nb_of_comp, nb_ve);
	heating_vis = alloc_2d_array<double>(nb_of_comp, nb_ve);
	heating_uv = alloc_2d_array<double>(nb_of_comp, nb_ve);
	heating_fuv = alloc_2d_array<double>(nb_of_comp, nb_ve);
	heating_cr = new double [nb_of_comp]; // no dependence on visual extinction;

	ISRF_Mathis1983 *is_field 
		= new ISRF_Mathis1983(); // without UV;

	ISRF_UV_Draine1978 *is_uv_field 
		= new ISRF_UV_Draine1978();

	CR_induced_UV_field *cr_field
		= new CR_induced_UV_field(cr_uv_flux);

	// the extinction law is set here, 31, 40, 55:
	fname = path + "dust/model_carb_silic_RV55_Draine2003.txt";
	dust_synthetic_Draine2003 *extinct_law 
		= new dust_synthetic_Draine2003(fname, verb = 0);

	// Note: the dust is irradiated from 4 pi here:
	a = 8*M_PI *PLANCK_CONSTANT *SPEED_OF_LIGHT *SPEED_OF_LIGHT;
	
	// radio  > 1 mm or < 10 cm-1
	// infrared 1 mm - 700 nm or 10 - 14290 cm-1; dust emission dominates at 0.6 mm - 5 um (Draine 2011)
	// visible 700 - 400 nm or 14290 - 25000 cm-1
	// ultraviolet 400 - 10 nm or 25000 - 
	for (j = 0; j < nb_ve; j++) 
	{
		// the H column density is calculated:
		f.l0 = extinction[j]*EXT_MAG_TO_OPT_DEPTH/extinct_law->opacity_perH(18182.); // opacity in cm2, 5500 A = 18182 cm-1;
		f.set(is_field, dust, extinct_law);
	
		for (i = 0; i < nb_of_comp; i++) 
		{
			f.i = i;
			heating_ra[i][j] = a *qromb<rf_func>(f, 0.1, 10., 1.e-5);
			heating_ir[i][j] = a *qromb<rf_func>(f, 10., 14290., 1.e-5);
			heating_vis[i][j] = a *qromb<rf_func>(f, 14290., 25000., 1.e-5);
			heating_uv[i][j] = a *qromb<rf_func>(f, 25000., 40816., 1.e-5);
		}

		f.set(is_uv_field, dust, extinct_law);
		for (i = 0; i < nb_of_comp; i++) 
		{
			f.i = i;		
			heating_fuv[i][j] = a *qromb<rf_func>(f, 40816., 109650., 1.e-5);	
		}
	}

	f.set(cr_field, dust, extinct_law);
	f.l0 = 0.; // no extinction for CR induced radiation;
	
	for (i = 0; i < nb_of_comp; i++) {
		f.i = i;	
		heating_cr[i] = a *qromb<rf_func>(f, 57143., 117650., 1.e-5); // 7.09-14.6 eV
	}

	delete is_field;
	delete is_uv_field;
	delete cr_field;
	delete extinct_law;
}

dust_heating_ISRF::~dust_heating_ISRF()
{
	delete [] extinction;
	delete [] heating_cr;
	free_2d_array(heating_ra);
	free_2d_array(heating_ir);
	free_2d_array(heating_vis);
	free_2d_array(heating_uv);
	free_2d_array(heating_fuv);
}

void dust_heating_ISRF::save_data(const std::string &path) const
{
	int i, j;
	double h;
	string fname;
	ofstream output;
	
	fname = path + "dust_isfr_heat";
	fname += name;
	fname += ".txt";
	output.open(fname.c_str(), ios_base::out);

	output << scientific;
	output.precision(3);
	output << left << "# The contrubition of interstellar infrared, optical, UV, CR induced UV radiation fields to dust heating, standard conditions" << endl 
		<< setw(13) << "# v_ext";
	
	for (i = 0; i < nb_of_comp; i++) {
		output << left << setw(13) << "total(erg/s)" << setw(13) << "ir(fraction) " << setw(13) << "vis" << setw(13) << "uv" 
			<< setw(13) << "fuv" << setw(13) << "cr" << "| ";
	}
	output << endl;

	for (j = 0; j < nb_ve; j++) 
	{
		output << left << setw(13) << extinction[j];
		for (i = 0; i < nb_of_comp; i++) 
		{
			h = heating_ra[i][j] + heating_ir[i][j] + heating_vis[i][j] + heating_uv[i][j] + heating_fuv[i][j] + heating_cr[i];
			output << left << setw(13) << h << setw(13) << heating_ir[i][j]/h << setw(13) << heating_vis[i][j]/h 
				<< setw(13) << heating_uv[i][j]/h << setw(13) << heating_fuv[i][j]/h << setw(13) << heating_cr[i]/h << "| ";
		}
		output << endl;
	}
	output.close();
}

void dust_heating_ISRF::get_heating(double vis_ext, double ir_field, double uv_field, double cr_ioniz, double *h, int nb) const
{
	if (nb != nb_of_comp) {
		cout << "Error in " << SOURCE_NAME << ": incorrect nb of dust components;";
		exit(1);
	}

	int j, l = 0, r = nb_ve-1;
	double h1, h2;
	while (r-l > 1)
	{
		j = l + ((r-l) >> 1);
		if (extinction[j] < vis_ext) 
			l = j;
		else r = j;
	}
	for (j = 0; j < nb_of_comp; j++) 
	{ 
		h1 = heating_ra[j][l] + ir_field*heating_ir[j][l] + heating_vis[j][l] + heating_uv[j][l] + uv_field*heating_fuv[j][l];
		h2 = heating_ra[j][l+1] + ir_field*heating_ir[j][l+1] + heating_vis[j][l+1] + heating_uv[j][l+1] + uv_field*heating_fuv[j][l+1];
		h[j] = h1 + (h2 - h1)/(extinction[l+1] - extinction[l]) *(vis_ext - extinction[l]) + cr_ioniz*heating_cr[j];
	}
}

double rf_func::operator() (double energy) const 
{
	// dust cross section for absorption, intensity normalized by 2hv/l^2;
	return energy*energy*energy *dust->components[i]->absorption(energy) 
		*rfield->get_intensity(energy) *exp(-l0*extinct_law->opacity_perH(energy)); 
}

//
// Functions
//

// For large carbonaceous and silicate grains:
// Note: to the output files one need to add the data: dust name, density and power of absorption dependence on wavelength;  
void reformat_dust_files(const string &path, string name1, string name2)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, nb_gr, nb_wl;
	double *wl_arr, *gr_arr, **abs_coeff, **scat_coeff, **cth; 

	string fname, str;
	ifstream input;
	ofstream output;

	fname = path + name1;
	input.open(fname.c_str());
	
	if (!input) {
		cout << " Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	for (i = 0; i < 3; i++) {
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	}
	
	input >> nb_gr;
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_wl;
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	gr_arr = new double [nb_gr];
	wl_arr = new double [nb_wl];

	abs_coeff = alloc_2d_array<double>(nb_gr, nb_wl);
	scat_coeff = alloc_2d_array<double>(nb_gr, nb_wl);
	cth = alloc_2d_array<double>(nb_gr, nb_wl);

	for (i = 0; i < nb_gr; i++)
	{
		input >> gr_arr[i];
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		for (j = 0; j < nb_wl; j++) 
		{
			input >> wl_arr[j];
			input >> abs_coeff[i][j] >> scat_coeff[i][j] >> cth[i][j];
		}
	}
	input.close();

	fname = path + name2;
	output.open(fname.c_str());
	output << scientific;
	output.precision(3);

	output << "# Here must be a reference" << endl << "# name, density (g/cm3), exponent for wavelength dependence at high wavelength;" << endl 
		<< "# grain radii is in cm; data: energy cm-1, abs, sca, cos theta" << endl 
		<< "name   denisty (g/cm3)   power" << endl;
	
	output << left << setw(5) << nb_wl << setw(5) << nb_gr << endl;
	output << left << setw(12) << " ";

	for (i = 0; i < nb_gr; i++) {
		// grain radius must be in cm:
		output << left << setw(36) << 1.e-4*gr_arr[i]; 
	}
	output << endl;
	
	for (j = 0; j < nb_wl; j++)
	{
		// conversion of wavelength from um to cm-1:
		output << left << setw(12) << 10000./wl_arr[j];
		for (i = 0; i < nb_gr; i++) {
			output << left << setw(12) << abs_coeff[i][j] << setw(12) << scat_coeff[i][j] << setw(12) << cth[i][j];
		}
		if (j < nb_wl-1) 
			output << endl;
	}
	output.close();

	delete [] gr_arr;
	delete [] wl_arr;
	free_2d_array(abs_coeff);
	free_2d_array(scat_coeff);
	free_2d_array(cth);
}

// For carbonaceous PAH grains, the data on PAH is supplemented by the data on carbonaceous grains at large grain radii > 100 A;
void reformat_dust_files_PAH(const string &path, string name1, string name2, string name3)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, k, nb1, nb2, nb_gr1, nb_gr2, nb_wl1, nb_wl2;
	
	double a;
	double *wl_arr1, *gr_arr1, *wl_arr2, *gr_arr2;
	double **abs_coeff1, **scat_coeff1, **cth1, **abs_coeff2, **scat_coeff2, **cth2; 

	string fname, str;
	ifstream input;
	ofstream output;

	// the first file is with PAH data:
	fname = path + name1;
	input.open(fname.c_str());

	if (!input) {
		cout << " Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	for (i = 0; i < 6; i++) {
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	}
	
	input >> nb_gr1;
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_wl1;
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	gr_arr1 = new double [nb_gr1];
	wl_arr1 = new double [nb_wl1];

	abs_coeff1 = alloc_2d_array<double>(nb_gr1, nb_wl1);
	scat_coeff1 = alloc_2d_array<double>(nb_gr1, nb_wl1);
	cth1 = alloc_2d_array<double>(nb_gr1, nb_wl1);

	for (i = 0; i < nb_gr1; i++)
	{
		input >> gr_arr1[i];
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		for (j = 0; j < nb_wl1; j++) 
		{
			input >> wl_arr1[j];
			input >> a; // extinction;
			input >> abs_coeff1[i][j] >> scat_coeff1[i][j] >> cth1[i][j];
		}
	}
	input.close();

	// the second file is with graphite data;
	fname = path + name2;
	input.open(fname.c_str());

	if (!input) {
		cout << " Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}

	for (i = 0; i < 3; i++) {
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	}
	
	input >> nb_gr2;
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_wl2;
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	gr_arr2 = new double [nb_gr2];
	wl_arr2 = new double [nb_wl2];

	abs_coeff2 = alloc_2d_array<double>(nb_gr2, nb_wl2);
	scat_coeff2 = alloc_2d_array<double>(nb_gr2, nb_wl2);
	cth2 = alloc_2d_array<double>(nb_gr2, nb_wl2);

	for (i = 0; i < nb_gr2; i++)
	{
		input >> gr_arr2[i];
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);

		for (j = 0; j < nb_wl2; j++) 
		{
			input >> wl_arr2[j];
			input >> abs_coeff2[i][j] >> scat_coeff2[i][j] >> cth2[i][j];
		}
	}

	for (nb1 = 0; (nb1 < nb_gr2) && (gr_arr1[nb_gr1-1] > 0.999*gr_arr2[nb1]); nb1++) 
	{;}
	// here grain radius is in microns, 1 micron = 1.e-4 cm; the maximal radius in question is 1.e-5 cm;
	for (nb2 = nb1; (nb2 < nb_gr2) && (gr_arr2[nb2] < 0.11); nb2++)
	{;}

	// the output data;
	fname = path + name3;
	output.open(fname.c_str());
	output << scientific;
	output.precision(3);

	output << "# Here must be a comment" << endl << "# name, density (g/cm3), exponent for wavelength dependence at high wavelength; " << endl 
		<< "# grain radii in cm; data: energy cm-1, abs, sca, cos theta; " << endl
		<<  "name   denisty (g/cm3)   power" << endl;
	output << left << setw(5) << nb_wl1 << setw(5) << nb_gr1 + nb2 - nb1 << endl;
	
	// grain radius must be in cm:
	output << left << setw(12) << " ";
	for (i = 0; i < nb_gr1; i++) {
		output << left << setw(36) << 1.e-4*gr_arr1[i]; 
	}
	for (i = nb1; i < nb2; i++) {
		output << left << setw(36) << 1.e-4*gr_arr2[i]; 
	}
	output << endl;
	
	for (j = 0; j < nb_wl1; j++)
	{
		// conversion of wavelength from um to cm-1:
		output << left << setw(12) << 10000./wl_arr1[j];
		// PAH data:
		for (i = 0; i < nb_gr1; i++) {
			output << left << setw(12) << abs_coeff1[i][j] << setw(12) << scat_coeff1[i][j] << setw(12) << cth1[i][j];
		}
		
		// wl_arr[] contains photon energy data;
		for (k = 1; (k < nb_wl2) && (wl_arr1[j] < wl_arr2[k]); k++) 
		{;}
		a = (wl_arr1[j] - wl_arr2[k-1])/(wl_arr2[k] -  wl_arr2[k-1]);

		// graphite data:
		for (i = nb1; i < nb2; i++) {
			output << left << setw(12) << abs_coeff2[i][k-1] + a*(abs_coeff2[i][k] - abs_coeff2[i][k-1])
				<< setw(12) << scat_coeff2[i][k-1] + a*(scat_coeff2[i][k] - scat_coeff2[i][k-1])
				<< setw(12) << cth2[i][k-1] + a*(cth2[i][k] - cth2[i][k-1]);
		}

		if (j < nb_wl1-1) 
			output << endl;
	}
	output.close();

	delete [] gr_arr1;
	delete [] gr_arr2;
	delete [] wl_arr1;
	delete [] wl_arr2;
	
	free_2d_array(abs_coeff1);
	free_2d_array(abs_coeff2);
	free_2d_array(scat_coeff1);
	free_2d_array(scat_coeff2);
	free_2d_array(cth1);
	free_2d_array(cth2);
}
