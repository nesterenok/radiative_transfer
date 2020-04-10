#pragma once
#include <math.h>
#include <string>

#define G_FUNC_PRECISION 1.e-6

// Hummer, Rybicki, ApJ 293, pp. 258-267, 1985;
// The functions are used to calculated the escape probability in the LVG (large velocity gradient) approximation; slab geometry;
// the result has to be multiplied by 0.5*g; 
class lvg_func_g {
public:
	double g, td; // td - dust optical depth from a point to the cloud boundary;
	double operator () (double mu) {
		return mu*mu*exp(-td/mu) *(1. - exp(-1./(g*mu*mu)));
	}
};

// Class containing the table with escape probability values;
// gamma = 1./( line absorp coeff * resonance length), delta = 1./( continuum absorp coeff * resonance length), 
// gamma > 0, delta > 0 for these data,
// across the resonance length the flow velocity shifts on the magnitude of the thermal velocity, 
// the function returns the value v, one-sided escape probability is p = 0.5(1-v), 	
class lvg_method_data {
protected:
	int		nb_g, nb_d;
	double	mu_c;
	double  *delta_arr, *gamma_arr, **p_arr;

public:
    // linear interpolation:
	virtual double get_esc_func(double gamma, double delta) const;
	
	// the function returns G(b, g, t) + G(b, g, T-t), see Hummer, Rybicki (1985);
	// d1 and d2 - optical depths in continuum to the cloud boundaries;
	double get_esc_func_g(double d1, double d2, double gamma) const;

	lvg_method_data(const std::string &path, std::string name, int verbosity = 1);
	virtual ~lvg_method_data();
};

class lvg_method_data_spline : public lvg_method_data 
{
protected:
    double ** p_d1_arr, ** p_d2_arr, ** p_d12_arr;

public:
    lvg_method_data_spline(const std::string& path, std::string name, int verbosity = 1);
    ~lvg_method_data_spline();

    // accurate interpolation:
    double get_esc_func(double gamma, double delta) const;
};

// not derived from previous classes,
class lvg_line_overlap_data
{
protected:
    int nb_d, nb_g, nb_gr, nb_dx;
    double *delta_arr, *gamma_arr, *gratio_arr, *dx_arr, ** p_arr;

public:
    // logarithmic interpolation for delta, linear interpolation for other parameters,
    // gamma, delta for first line, delta = d1, gamma = g1, gamma_ratio = g2/g1, dx = (fr1 - fr2) / dfr, linear interpolation,
    // gamma > 0, delta > 0 (for gamma < 0, delta < 0, one has to  gamma -> -gamma, delta -> -delta, dx -> -dx)
    double get_esc_func(double gamma, double delta, double gamma_ratio, double delta_x) const;
    double get_dx_lim() const { return fabs(dx_arr[0]); } // the array values may not be symmetric,

    lvg_line_overlap_data(const std::string& path, std::string name, int verbosity = 1);
    ~lvg_line_overlap_data();
};

/*
// Sobovev, Soviet Astronomy 1, p.678, 1957; Hummer, Rybicki, ApJ 254, p. 767-779, 1982;
// the function has to be integrated on mu to derive one-sided photon escape probability from a slab;
// large velocity gradient approximation, no continuum opacity; the result has to be multiplied by 0.5*g;
class lvg_func_slab {
public:
	double g;
	double operator() (double mu) {
		return mu * mu * (1. - exp(-1. / (g * mu * mu)));
	}
};
*/
