
//
// 07.03.2017. Check for errors.
// 31.08.2017.  Check for errors.
// 30.01.2018. Check for errors.

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <fstream>
#include <iostream>

#include "constants.h"
#include "utils.h"
#include "integration.h"
#include "photoelectric_emission.h"

#define MAX_TEXT_LINE_WIDTH 240
#define SOURCE_NAME "photoelectric_emission.cpp"
using namespace std;

//
// Base class for calculation of photoelectric rates
//

photoel_emission::photoel_emission(int v) : rfield(0), verbosity(v) {
	length_el = 1.e-7; // // electron escape length in cm;
}

photoel_emission::~photoel_emission() 
{
	delete [] ph_en_arr;
	delete [] ph_en_arr2;
	delete [] grain_rad_arr;
	delete [] abs_length;
	free_2d_array(abs_coeff);
}

double photoel_emission::get_abs_length(double eph) const
{
	if (eph < ph_en_arr2[0])
		return abs_length[0];

	else if (eph > ph_en_arr2[nb_ph_en2-1])
		return abs_length[nb_ph_en2-1];
	
	int k;
	locate_index(ph_en_arr2, nb_ph_en2, eph, k);
	return abs_length[k] + (abs_length[k+1] - abs_length[k])*(eph - ph_en_arr2[k])/(ph_en_arr2[k+1] - ph_en_arr2[k]);
}

void photoel_emission::set_radius(double r)
{
	radius = r;
	locate_index(grain_rad_arr, nb_gr_rad, radius, index_r);
	
	if (index_r < 0) 
	{
		index_r = 0;
		t = 0.;
	}
	else if (index_r >= nb_gr_rad-1) 
	{
		index_r = nb_gr_rad-2;
		t = 1.;
	}
	else t = (radius - grain_rad_arr[index_r])/(grain_rad_arr[index_r+1] - grain_rad_arr[index_r]);
}

double photoel_emission::get_abs_cs(double eph) const
{
	int l;
	double u;

	locate_index(ph_en_arr, nb_ph_en, eph, l);
	if (l < 0) {
		l = 0;
		u = 0.;
	}
	else if (l >= nb_ph_en-1) {
		l = nb_ph_en-2;
		u = 1.;
	}
	else u = (eph - ph_en_arr[l])/(ph_en_arr[l+1] - ph_en_arr[l]);

	return abs_coeff[l][index_r]*(1.-u)*(1.-t) + abs_coeff[l+1][index_r]*u*(1.-t) + abs_coeff[l][index_r+1]*(1.-u)*t + abs_coeff[l+1][index_r+1]*u*t;
}

//
// Calculation of photoelectric yields for carbonaceous dust
//

void photoel_emission_carbon_dust::get_limit_charges(double radius, int & zmin, int & zmax) const
{
	// different for carbonaceous and silicate grains:
	zmin = (int) ((-3.9e+8*radius - 0.12e+16*radius*radius - 2.)/14.4);

	// common for both grain types:
	zmax = (int) ( ((13.6*EV_TO_CM_INVERSE - work_f)/(14.4*EV_TO_CM_INVERSE) *1.e+8*radius + 0.5 - 0.3e-8/radius)
				/(1 + 0.3e-8/radius) );
}

photoel_emission_carbon_dust::photoel_emission_carbon_dust(const string &path, int verb) : photoel_emission(verb)
{
	char text_line[MAX_TEXT_LINE_WIDTH];	
	int i, j;
	double a;

	string fname, str;
	ifstream input;

	// graphite work function:
	work_f = 4.4*EV_TO_CM_INVERSE;

	// data for PAH
	fname = path + "dust/dust_PAHneut_Draine2001.txt";
	input.open(fname.c_str(), ios_base::in);
	
	if (!input) {
		cout << " Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> str >> a >> a;
	input >> nb_ph_en3 >> nb_pah_rad;

	ph_en_arr3 = new double [nb_ph_en3];
	pah_rad_arr = new double [nb_pah_rad];
	abs_coeff_p = alloc_2d_array<double>(nb_ph_en3, nb_pah_rad);

	for (i = 0; i < nb_pah_rad; i++) {
		input >> pah_rad_arr[i];
	}
	i = 0;
	while (!input.eof() && i < nb_ph_en3)
	{
		input >> ph_en_arr3[i]; // in cm-1
		for (j = 0; j < nb_pah_rad; j++) 
		{
			input >> abs_coeff_p[i][j] >> a >> a;
			abs_coeff_p[i][j] *= M_PI*pah_rad_arr[j]*pah_rad_arr[j];
		}
		i++;
	}
	input.close();

	// data for large carbonaceois grains, data in the file is for grain sizes > 10 A = 1.e-7 cm;
	fname = path + "dust/dust_graphite_Draine1993.txt";
	input.open(fname.c_str(), ios_base::in);
	
	if (!input) {
		cout << " Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> str >> a >> a;
	input >> nb_ph_en >> nb_gr_rad;

	ph_en_arr = new double [nb_ph_en];
	grain_rad_arr = new double [nb_gr_rad];
	abs_coeff = alloc_2d_array<double>(nb_ph_en, nb_gr_rad);

	for (i = 0; i < nb_gr_rad; i++) {
		input >> grain_rad_arr[i];
	}
	i = 0;
	while (!input.eof() && i < nb_ph_en)
	{
		input >> ph_en_arr[i]; // in cm-1
		for (j = 0; j < nb_gr_rad; j++) 
		{
			input >> abs_coeff[i][j] >> a >> a;
			abs_coeff[i][j] *= M_PI*grain_rad_arr[j]*grain_rad_arr[j];
		}
		i++;
	}
	input.close();

	fname = path + "dust/abs_length_carbonaceous.txt";
	input.open(fname.c_str(), ios_base::in);
	
	if (!input) {
		cout << " Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_ph_en2;
	ph_en_arr2 = new double [nb_ph_en2];
	abs_length = new double [nb_ph_en2];

	for (i = 0; i < nb_ph_en2; i++) 
	{
		// in file photon wavelength is in um, absorption length in cm;
		input >> a >> abs_length[i];  
		// photon energy must be in cm-1; the photon energy is arranged in ascending order;
		ph_en_arr2[i] = 1.e+4/a;
	}
	input.close();
}

photoel_emission_carbon_dust::~photoel_emission_carbon_dust()
{
	delete [] ph_en_arr3;
	delete [] pah_rad_arr;
	free_2d_array(abs_coeff_p);
}

void photoel_emission_carbon_dust::set_radius(double r)
{
	if (r >= 1.e-6) { // in cm, 1.e-6 cm = 100 A
		// is not PAH;
		photoel_emission::set_radius(r);
	}
	else
	{
		// is PAH;
		radius = r;
		locate_index(pah_rad_arr, nb_pah_rad, radius, index_r);
	
		if (index_r < 0) 
		{
			index_r = 0;
			t = 0.;
		}
		else if (index_r >= nb_pah_rad-1) 
		{
			index_r = nb_pah_rad-2;
			t = 1.;
		}
		else t = (radius - pah_rad_arr[index_r])/(pah_rad_arr[index_r+1] - pah_rad_arr[index_r]);
	}
}

double photoel_emission_carbon_dust::get_abs_cs(double eph) const
{
	if (radius >= 1.e-6) { // in cm
		return photoel_emission::get_abs_cs(eph);
	}
	else
	{
		int l;
		double u;

		locate_index(ph_en_arr3, nb_ph_en3, eph, l);
		if (l < 0) {
			l = 0;
			u = 0.;
		}
		else if (l >= nb_ph_en3-1) {
			l = nb_ph_en3-2;
			u = 1.;
		}
		else u = (eph - ph_en_arr3[l])/(ph_en_arr3[l+1] - ph_en_arr3[l]);

		return abs_coeff_p[l][index_r]*(1.-u)*(1.-t) + abs_coeff_p[l+1][index_r]*u*(1.-t) 
			+ abs_coeff_p[l][index_r+1]*(1.-u)*t + abs_coeff_p[l+1][index_r+1]*u*t;
	}
}

double photoel_emission_carbon_dust::get_yield(int z, double radius, double eph) const
{
	double a, b, y0, y1, y2, yield, eph_thr, emin, elow, ehigh, theta;
	
	if (z < 0)
		emin = -(z + 1)*ELECTRON_CHARGE*ELECTRON_CHARGE /(radius*(pow(2.7e-7/radius, 0.75) + 1.)) *ERG_TO_CM_INVERSE;
	else emin = 0.;

	eph_thr = work_f + ELECTRON_CHARGE*ELECTRON_CHARGE/radius*(z + 0.5 + 3.e-9*(z + 2)/radius)*ERG_TO_CM_INVERSE;
	if (z < -1)
		eph_thr += emin;

	if (eph_thr < 0.) {
		eph_thr = 0.;
	}
	
	if (eph < eph_thr)
		return 0.;

	if (z >= 0)
		elow = -(z + 1)*ELECTRON_CHARGE*ELECTRON_CHARGE/radius *ERG_TO_CM_INVERSE;
	else elow = emin;

	if (z >= 0) 
	{
		theta = eph - eph_thr - elow;
		ehigh = eph - eph_thr;
		y2 = ehigh*ehigh*(ehigh - 3.*elow)/((ehigh - elow)*(ehigh - elow)*(ehigh - elow));
	}
	else {
		theta = eph - eph_thr;	
		ehigh = emin + eph - eph_thr;
		y2 = 1.;
	}

	a = pow(theta/work_f, 5.);
	y0 = 9.e-3*a/( 1. + 3.7e-2*a );

	b = radius/get_abs_length(eph);
	a = b + radius/length_el;
		
	y1 = b*b*(a*a - 2.*a + 2. - 2.*exp(-a)) /( a*a*(b*b - 2.*b + 2. - 2.*exp(-b)) );

	yield = y2;
	if (y0*y1 < 1.)
		yield *= y0*y1;

	return yield;
}

void photoel_emission_carbon_dust::get_phel_emiss(int z, double radius, double &phem_rate, double &av_energy)
{
	int i, nb;
	double a, b, c, de, ea, theta, emin, elow, ehigh, eph, eph_thr, y0, y1, y2, yield;

	set_radius(radius);
	
	if (z < 0)
		emin = -(z + 1)*ELECTRON_CHARGE*ELECTRON_CHARGE /(radius*(pow(2.7e-7/radius, 0.75) + 1.)) *ERG_TO_CM_INVERSE;
	else emin = 0.;

	eph_thr = work_f + ELECTRON_CHARGE*ELECTRON_CHARGE/radius*(z + 0.5 + 3.e-9*(z + 2)/radius)*ERG_TO_CM_INVERSE;
	if (z < -1)
		eph_thr += emin;

	if (eph_thr < 0.) {
		eph_thr = 0.;
		if (verbosity)
			cout << "Warning: for valence band electron threshold energy is < 0." << endl;
	}
	
	if (z >= 0)
		elow = -(z + 1)*ELECTRON_CHARGE*ELECTRON_CHARGE/radius *ERG_TO_CM_INVERSE;
	else elow = emin;

	// Valence band electrons;
	phem_rate = av_energy = 0.;
	if (eph_thr < 13.6*EV_TO_CM_INVERSE)
	{
		// trapezoidal method is employed for integration by energy, with the fixed number of intervals; 
		nb = 1000;
		de = (13.6*EV_TO_CM_INVERSE - eph_thr)/nb;	
		eph = eph_thr;

		for (i = 0; i <= nb; i++)
		{
			if (z >= 0) 
			{
				theta = eph - eph_thr - elow;
				ehigh = eph - eph_thr;
				y2 = ehigh*ehigh*(ehigh - 3.*elow)/((ehigh - elow)*(ehigh - elow)*(ehigh - elow));
			}
			else {
				theta = eph - eph_thr;	
				ehigh = emin + eph - eph_thr;
				y2 = 1.;
			}

			a = pow(theta/work_f, 5.);
			y0 = 9.e-3*a/( 1. + 3.7e-2*a );

			b = radius/get_abs_length(eph);
			a = b + radius/length_el;
		
			y1 = b*b*(a*a - 2.*a + 2. - 2.*exp(-a)) /( a*a*(b*b - 2.*b + 2. - 2.*exp(-b)) );

			yield = y2;
			if (y0*y1 < 1.)
				yield *= y0*y1;

			// the intensity provided by function must be multiplied by 2h*v/lambda^2 *4pi and divided by hv, 
			// integrated by energy in cm-1 - must be multiplied by speed of light; see below the normalization; 
			// the rate of electron production:
			if (eph > DBL_EPSILON)
				c = yield *get_abs_cs(eph) *rfield->get_intensity(eph) *eph*eph;
			else c = 0.;

			if (i == 0 || i == nb)
				c *= 0.5;

			phem_rate += c;
		
			// the calculation of electron energy based on distribution by Weingartner & Draine (2001),
			// the energy of the electron with respect to infinity;
			if (z >= 0) {
				av_energy += 0.5*ehigh*(ehigh - 2.*elow)/(ehigh - 3.*elow) *c;
			}
			else { 
				if (fabs(ehigh - elow) > DBL_EPSILON) // ?
					av_energy += 0.5*(ehigh*ehigh*ehigh*(ehigh - 2.*elow) - elow*elow*elow*(elow - 2.*ehigh))
						/((ehigh - elow)*(ehigh - elow)*(ehigh - elow)) *c;
				else 
					av_energy += elow *c;
			}
			eph += de;
		}
		phem_rate *= de;
		av_energy *= de;
	}
	
	// Photodetachment electrons;
	y0 = y1 = 0.; // auxulary variables;
	if (z < 0) 
	{
		// EA(Z+1)
		ea = work_f + ((z + 0.5) - 4.e-8/(radius + 7.e-8)) *ELECTRON_CHARGE*ELECTRON_CHARGE/radius *ERG_TO_CM_INVERSE;
		eph_thr = ea + emin;

		if (eph_thr < 0.) {
			eph_thr = 0.;
			if (verbosity)
				cout << "Warning: photodetachment threshold energy is < 0." << endl;
		}
		
		if (eph_thr < 13.6*EV_TO_CM_INVERSE)
		{
			de = (13.6*EV_TO_CM_INVERSE - eph_thr)/nb;	
			eph = eph_thr;

			for (i = 0; i <= nb; i++)
			{
				a = (eph - eph_thr)/(3.*EV_TO_CM_INVERSE);
				b = 1. + a*a/3.;
				
				if (eph > DBL_EPSILON)
					c = 1.2e-17*abs(z)*a/(b*b) *rfield->get_intensity(eph) *eph*eph;
				else c = 0.;

				if (i == 0 || i == nb)
					c *= 0.5;

				y0 += c;
				y1 += c*(eph - eph_thr + emin);
				eph += de;
			}
			y0 *= de;
			y1 *= de;
		}
	}
	phem_rate += y0;
	av_energy += y1;

	// speed of light constant appears from the change of integration variable Hz -> cm-1:
	phem_rate *= 8.*M_PI*SPEED_OF_LIGHT;
	av_energy *= 8.*M_PI*SPEED_OF_LIGHT;

	if (phem_rate < MIN_PHOTOEMISSION_RATE) // some arbitrary small value;
		phem_rate = av_energy = 0.;
	else av_energy /= phem_rate;
}

//
// Calculation of photoelectric yields for silicate dust
//

void photoel_emission_silicate_dust::get_limit_charges(double radius, int & zmin, int & zmax) const
{
	// different from carbonaceous grains:
	zmin = (int) ((-2.5e+8*radius - 0.07e+16*radius*radius - 8.)/14.4);

	// common for both grain types:
	zmax = (int) ( ((13.6*EV_TO_CM_INVERSE - work_f)/(14.4*EV_TO_CM_INVERSE) *1.e+8*radius + 0.5 - 0.3e-8/radius)
				/(1 + 0.3e-8/radius) );
}

photoel_emission_silicate_dust::photoel_emission_silicate_dust(const string &path, int verb) : photoel_emission(verb)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j;
	double a;

	string fname, str;
	ifstream input;

	// silicate work function:
	work_f = 8.*EV_TO_CM_INVERSE;

	fname = path + "dust/dust_silicate_Draine2001.txt";
	input.open(fname.c_str(), ios_base::in);
	
	if (!input) {
		cout << " Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> str >> a >> a;
	input >> nb_ph_en >> nb_gr_rad;

	ph_en_arr = new double [nb_ph_en];
	grain_rad_arr = new double [nb_gr_rad];
	abs_coeff = alloc_2d_array<double>(nb_ph_en, nb_gr_rad);

	for (i = 0; i < nb_gr_rad; i++) {
		input >> grain_rad_arr[i];
	}
	
	i = 0;
	while (!input.eof() && i < nb_ph_en)
	{
		input >> ph_en_arr[i]; // in cm-1
		for (j = 0; j < nb_gr_rad; j++) 
		{
			input >> abs_coeff[i][j] >> a >> a;
			abs_coeff[i][j] *= M_PI*grain_rad_arr[j]*grain_rad_arr[j];
		}
		i++;
	}
	input.close();

	fname = path + "dust/abs_length_silicate.txt";
	input.open(fname.c_str(), ios_base::in);
	
	if (!input) {
		cout << " Error in " << SOURCE_NAME << ": can't open " << fname << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> nb_ph_en2;
	ph_en_arr2 = new double [nb_ph_en2];
	abs_length = new double [nb_ph_en2];

	for (i = 0; i < nb_ph_en2; i++) 
	{
		input >> a >> abs_length[i]; // length in cm;
		// photon energy must be in cm-1, the photon energy is arranged in ascending order;
		ph_en_arr2[i] = 1.e+4/a; 
	}
	input.close();
}

double photoel_emission_silicate_dust::get_yield(int z, double radius, double eph) const
{
	double a, b, y0, y1, y2, yield, theta, eph_thr, emin, elow, ehigh;
	if (z < 0)
		emin = -(z + 1)*ELECTRON_CHARGE*ELECTRON_CHARGE /(radius*(pow(2.7e-7/radius, 0.75) + 1.)) *ERG_TO_CM_INVERSE;
	else emin = 0.;

	eph_thr = work_f + ELECTRON_CHARGE*ELECTRON_CHARGE/radius *(z + 0.5 + 3.e-9*(z + 2)/radius)*ERG_TO_CM_INVERSE;
	if (z < -1)
		eph_thr += emin;

	if (eph_thr < 0.) {
		eph_thr = 0.;
	}

	if (eph < eph_thr)
		return 0.;

	if (z >= 0)
		elow = -(z + 1)*ELECTRON_CHARGE*ELECTRON_CHARGE/radius *ERG_TO_CM_INVERSE;
	else elow = emin;

	if (z >= 0) 
	{
		theta = eph - eph_thr - elow;		
		ehigh = eph - eph_thr;
		y2 = ehigh*ehigh*(ehigh - 3.*elow)/((ehigh - elow)*(ehigh - elow)*(ehigh - elow));
	}
	else {
		theta = eph - eph_thr;
		ehigh = emin + eph - eph_thr;
		y2 = 1.;
	}

	a = theta/work_f;
	y0 = 0.5*a/( 1. + 5.*a );

	b = radius/get_abs_length(eph);
	a = b + radius/length_el;
		
	y1 = b*b*(a*a - 2.*a + 2. - 2.*exp(-a)) /( a*a*(b*b - 2.*b + 2. - 2.*exp(-b)) );

	yield = y2;
	if (y0*y1 < 1.)
		yield *= y0*y1;
	
	return yield;
}

void photoel_emission_silicate_dust::get_phel_emiss(int z, double radius, double &phem_rate, double &av_energy)
{
	int i, nb;
	double a, b, c, de, ea, theta, emin, elow, ehigh, eph, eph_thr, y0, y1, y2, yield;

	set_radius(radius);
	
	if (z < 0)
		emin = -(z + 1)*ELECTRON_CHARGE*ELECTRON_CHARGE /(radius*(pow(2.7e-7/radius, 0.75) + 1.)) *ERG_TO_CM_INVERSE;
	else emin = 0.;

	eph_thr = work_f + ELECTRON_CHARGE*ELECTRON_CHARGE/radius *(z + 0.5 + 3.e-9*(z + 2)/radius)*ERG_TO_CM_INVERSE;
	if (z < -1)
		eph_thr += emin;

	if (eph_thr < 0.) {
		eph_thr = 0.;
		if (verbosity)
			cout << "Warning: for valence band electron threshold energy is < 0." << endl;
	}

	if (z >= 0)
		elow = -(z + 1)*ELECTRON_CHARGE*ELECTRON_CHARGE/radius *ERG_TO_CM_INVERSE;
	else elow = emin;

	// Valence band electrons;
	phem_rate = av_energy = 0.;
	if (eph_thr < 13.6*EV_TO_CM_INVERSE)
	{
		// trapezoidal method is employed for integration by energy, with the fixed number of intervals; 
		nb = 1000;
		de = (13.6*EV_TO_CM_INVERSE - eph_thr)/nb;
		eph = eph_thr;

		for (i = 0; i <= nb; i++)
		{
			if (z >= 0) 
			{
				theta = eph - eph_thr - elow;		
				ehigh = eph - eph_thr;
				y2 = ehigh*ehigh*(ehigh - 3.*elow)/((ehigh - elow)*(ehigh - elow)*(ehigh - elow));
			}
			else {
				theta = eph - eph_thr;
				ehigh = emin + eph - eph_thr;
				y2 = 1.;
			}

			a = theta/work_f;
			y0 = 0.5*a/( 1. + 5.*a );

			b = radius/get_abs_length(eph);
			a = b + radius/length_el;
		
			y1 = b*b*(a*a - 2.*a + 2. - 2.*exp(-a)) /( a*a*(b*b - 2.*b + 2. - 2.*exp(-b)) );

			yield = y2;
			if (y0*y1 < 1.)
				yield *= y0*y1;

			// the intensity provided by function must be multiplied by 2h*v/lambda^2 *4pi and divided by hv, 
			// integrated by energy in cm-1 - must be multiplied by speed of light; see below the normalization; 
			// the rate of electron production:
			if (eph > DBL_EPSILON)
				c = yield *get_abs_cs(eph) *rfield->get_intensity(eph) *eph*eph;
			else c = 0.;

			if (i == 0 || i == nb)
				c *= 0.5;
		
			phem_rate += c;
			if (z >= 0) {
				av_energy += 0.5*ehigh*(ehigh - 2.*elow)/(ehigh - 3.*elow) *c;
			}
			else {
				if (fabs(ehigh - elow) > DBL_EPSILON)
					av_energy += 0.5*(ehigh*ehigh*ehigh*(ehigh - 2.*elow) - elow*elow*elow*(elow - 2.*ehigh))
						/((ehigh - elow)*(ehigh - elow)*(ehigh - elow)) *c;
				else av_energy += elow *c;
			}
			eph += de;
		}
		phem_rate *= de;
		av_energy *= de;
	}

	// Photodetachment electrons;
	y0 = y1 = 0.; // auxulary variables;
	if (z < 0) 
	{
		// EA(Z+1)
		ea = work_f - 5.*EV_TO_CM_INVERSE + (z + 0.5)*ELECTRON_CHARGE*ELECTRON_CHARGE/radius *ERG_TO_CM_INVERSE;
		eph_thr = ea + emin;

		if (eph_thr < 0.) {
			eph_thr = 0.;
			if (verbosity)
				cout << "Warning: photodetachment threshold energy is < 0." << endl;
		}

		if (eph_thr < 13.6*EV_TO_CM_INVERSE)
		{
			de = (13.6*EV_TO_CM_INVERSE - eph_thr)/nb;	
			eph = eph_thr;

			for (i = 0; i <= nb; i++)
			{
				a = (eph - eph_thr)/(3.*EV_TO_CM_INVERSE);
				b = 1. + a*a/3.;
				
				if (eph > DBL_EPSILON)
					c = 1.2e-17*abs(z)*a/(b*b) *rfield->get_intensity(eph) *eph*eph;
				else c = 0.;

				if (i == 0 || i == nb)
					c *= 0.5;

				y0 += c;
				y1 += c*(eph - eph_thr + emin);
				eph += de;
			}
			y0 *= de;
			y1 *= de;
		}
	}
	phem_rate += y0;
	av_energy += y1;

	// speed of light constant appears from the change of integration variable Hz -> cm-1:
	phem_rate *= 8.*M_PI*SPEED_OF_LIGHT;
	av_energy *= 8.*M_PI*SPEED_OF_LIGHT;

	if (phem_rate < MIN_PHOTOEMISSION_RATE) // some arbitrary small value;
		phem_rate = av_energy = 0.;
	else av_energy /= phem_rate;
}

// Only UV fields are considered:
void calc_grain_photoelectron_rates(const std::string &input_path, double cr_uv_flux)
{
	const int z_lim_low = -100;
	const int z_lim_up = 100;

	int i, j, z, zmin, zmax, nb_r, nb_eph;
	double radius, rmin, rmax, sc, rate, energy;

	string fname, output_path = "";
	ofstream output;

	// nb of energies in the calculation of photoelectric yields:
	nb_eph = 30;

	// minimal and maximal radii in cm;
	rmin = 3.5e-8; // 1 A = 1.e-8 cm
	rmax = 1.e-4; // 1 um = 1.e-4 cm
	
	sc = pow(10, 0.25);
	nb_r = (int) (log(rmax/rmin)/log(sc)) + 2;

	photoel_emission_carbon_dust *phem_carb =
		new photoel_emission_carbon_dust(input_path);

	photoel_emission_silicate_dust *phem_sil =
		new photoel_emission_silicate_dust(input_path);

	ISRF_UV_Draine1978 *rf1 
		= new ISRF_UV_Draine1978();

	CR_induced_UV_field *rf2 
		= new CR_induced_UV_field(cr_uv_flux);
	
	fname = output_path + "dust_photoelectric_rates.txt";
	output.open(fname.c_str());

	output << scientific;
	output.precision(2);

	output << "# Photoelectric yields for neutral carbonaceous grains:" << endl
		<< "on x axis - radius in cm, on y axis - energy in cm-1" << endl;

	output << left << setw(11) << "";
	for (radius = rmin, i = 0; i < nb_r; radius *= sc, i++) {
		output << left << setw(11) << radius;
	}
	output << endl;
	
	for (j = 0; j < nb_eph; j++)
	{
		energy = (13.6*EV_TO_CM_INVERSE*j)/nb_eph;
		output << left << setw(11) << energy;

		for (radius = rmin, i = 0; i < nb_r; radius *= sc, i++) {
			output << left << setw(11) << phem_carb->get_yield(z = 0, radius, energy);
		}
		output << endl;
	}
	output << endl;

	output << "# Photoelectric yields for neutral silicate grains:" << endl
		<< "on x axis - radius in cm, on y axis - energy in cm-1" << endl;

	output << left << setw(11) << "";
	for (radius = rmin, i = 0; i < nb_r; radius *= sc, i++) {
		output << left << setw(11) << radius;
	}
	output << endl;
	
	for (j = 0; j < nb_eph; j++)
	{
		energy = (13.6*EV_TO_CM_INVERSE*j)/nb_eph;
		output << left << setw(11) << energy;

		for (radius = rmin, i = 0; i < nb_r; radius *= sc, i++) {
			output << left << setw(11) << phem_sil->get_yield(z = 0, radius, energy);
		}
		output << endl;
	}
	output << endl;

	output << "# Photoelectric rates on grains;" << endl
		<< "# on horizantol axis - charge, on vertical axis - grain radius (in cm); values - rate (s-1) and electron energy (erg);" << endl 
		<< "# first table - photoeffect induced by interstellar UV radiation field on carbonaceous grains;" << endl;

	output << left << setw(5) << nb_r << setw(5) << z_lim_up - z_lim_low + 1 << endl;

	output << left << setw(11) << "";
	for (z = z_lim_low; z <= z_lim_up; z++) {
		output << left << setw(11) << z << setw(11) << z;
	}
	output << endl;

	phem_carb->set_rfield(rf1);
	
	for (radius = rmin, i = 0; i < nb_r; radius *= sc, i++) 
	{
		phem_carb->get_limit_charges(radius, zmin, zmax);
		if (zmin < z_lim_low)
			zmin = z_lim_low;
			
		if (zmax > z_lim_up)
			zmax = z_lim_up;

		output << left << setw(11) << radius;
		for (z = z_lim_low; z <= z_lim_up; z++) 
		{
			if (z >= zmin && z <= zmax) 
			{
				phem_carb->get_phel_emiss(z, radius, rate, energy);
				output << left << setw(11) << rate << setw(11) << energy*CM_INVERSE_TO_ERG;
			}
			else {
				output << left << setw(11) << "0." << setw(11) << "0.";
			}
		}
		output << endl;
	}
	
	phem_carb->set_rfield(rf2);
	output << endl << "# Second table - photoeffect induced by CR induced radiation field on cabonaceous grains;" << endl;

	for (radius = rmin, i = 0; i < nb_r; radius *= sc, i++) 
	{
		phem_carb->get_limit_charges(radius, zmin, zmax);

		if (zmin < z_lim_low)
			zmin = z_lim_low;
			
		if (zmax > z_lim_up)
			zmax = z_lim_up;

		output << left << setw(11) << radius;
		for (z = z_lim_low; z <= z_lim_up; z++) 
		{
			if (z >= zmin && z <= zmax) 
			{			
				phem_carb->get_phel_emiss(z, radius, rate, energy);
				output << left << setw(11) << rate << setw(11) << energy*CM_INVERSE_TO_ERG;
			}
			else {
				output << left << setw(11) << "0." << setw(11) << "0.";
			}
		}
		output << endl;
	}

	phem_sil->set_rfield(rf1);
	output << endl << "# Third table - photoeffect induced by interstellar UV radiation field on silicate grains;" << endl
		<< "# for silicate grains, the data on absorption coeffisient exist for radii 10 A <= r <= 1e+5 A = 10 um" << endl;
	
	for (radius = rmin, i = 0; i < nb_r; radius *= sc, i++) 
	{
		phem_sil->get_limit_charges(radius, zmin, zmax);

		if (zmin < z_lim_low)
			zmin = z_lim_low;
			
		if (zmax > z_lim_up)
			zmax = z_lim_up;

		output << left << setw(11) << radius;
		for (z = z_lim_low; z <= z_lim_up; z++) 
		{
			if (z >= zmin && z <= zmax) 
			{
				phem_sil->get_phel_emiss(z, radius, rate, energy);
				output << left << setw(11) << rate << setw(11) << energy*CM_INVERSE_TO_ERG;
			}
			else {
				output << left << setw(11) << "0." << setw(11) << "0.";
			}
		}
		output << endl;
	}
	
	phem_sil->set_rfield(rf2);
	output << endl << "# Fourth table - photoeffect induced by CR induced radiation field on silicate grains;" << endl;

	for (radius = rmin, i = 0; i < nb_r; radius *= sc, i++) 
	{
		phem_sil->get_limit_charges(radius, zmin, zmax);

		if (zmin < z_lim_low)
			zmin = z_lim_low;
			
		if (zmax > z_lim_up)
			zmax = z_lim_up;

		output << left << setw(11) << radius;
		for (z = z_lim_low; z <= z_lim_up; z++) 
		{
			if (z >= zmin && z <= zmax) 
			{
				phem_sil->get_phel_emiss(z, radius, rate, energy);
				output << left << setw(11) << rate << setw(11) << energy*CM_INVERSE_TO_ERG;
			}
			else {
				output << left << setw(11) << "0." << setw(11) << "0.";
			}
		}
		output << endl;
	}
	output.close();

	delete rf1;
	delete rf2;
	delete phem_carb;
	delete phem_sil;
}
