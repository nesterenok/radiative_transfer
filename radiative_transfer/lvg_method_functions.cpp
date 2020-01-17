
//
// 05.09.2017. Check for errors;
// 13.01.2020. Check for errors, line overlap was added;

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cmath>

#include "lvg_method_functions.h"
#include "integration.h"
#include "interpolation.h"
#include "special_functions.h"

#define MAX_TEXT_LINE_WIDTH 240
using namespace std;

lvg_method_data::lvg_method_data(const string &path, string name, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i,j;
	
	string file_name;
	ifstream input;

    file_name = path + name;
	input.open(file_name.c_str(), ios_base::in);

	if (!input.is_open()) {
		cout << "Error in lvg_method_functions.cpp: can't open " << file_name << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> nb_d >> nb_g;

	// Delta is the 1st variable, gamma - the 2nd;
	delta_arr = new double [nb_d];
	gamma_arr = new double [nb_g];
	
	p_arr = alloc_2d_array<double>(nb_d, nb_g);
	for (i = 0; i < nb_g; i++) {
		input >> gamma_arr[i];
	}

	for (i = 0; i < nb_d; i++) {
		input >> j >> delta_arr[i];
		for (j = 0; j < nb_g; j++) {
			input >> p_arr[i][j];
		}
	}
	input.close();

	// the minimal cosine of the angle used in integration;
	mu_c = 1.e-6;

	if (verbosity)
		cout << "The data on LVG method have been initialized." << endl;
}

lvg_method_data::~lvg_method_data()
{
	delete [] gamma_arr;
	delete [] delta_arr;
	free_2d_array(p_arr);
}

// linear interpolation for the escape probability function p(d,g);
double lvg_method_data::get_esc_func(double gamma, double delta) const
{
	int k, l;
	double t, u, escf(0.);
	
	locate_index(delta_arr, nb_d, delta, k);
	locate_index(gamma_arr, nb_g, gamma, l);

	if (k < 0) {
		t = 0.;
		k = 0;
	}
	else if (k > nb_d-2) {
		t = 1.;
		k = nb_d-2;
	}
	else t = (delta - delta_arr[k])/(delta_arr[k+1] - delta_arr[k]);
	
	if (l < 0) { 
		l = 0;
		u = 0.;
	}
	else if (l > nb_g-2) { 
		l = nb_g-2;
		u = 1.;
	}
	else u = (gamma - gamma_arr[l])/(gamma_arr[l+1] - gamma_arr[l]);
	
	escf = p_arr[k][l]*(1.-t)*(1.-u) + p_arr[k+1][l]*t*(1.-u) + p_arr[k][l+1]*(1.-t)*u + p_arr[k+1][l+1]*u*t;

	if (escf > 1.) 
		escf = 1.;
	else if (escf < 0.) 
		escf = 0.;

	return escf;
}

double lvg_method_data::get_esc_func_g(double d1, double d2, double gamma) const
{
	double escf;	
	lvg_func_g lvg_fc;

	if (gamma < 0.0001) {
		escf = 0.5*gamma *(expint(4, d1) + expint(4, d2));
	}
	else
	{
		lvg_fc.g = gamma;
	// the first condition - the gamma must be sufficiently large; the second - the convergence criterion on series;
		if (mu_c*mu_c*gamma > 1. && d1*d1*gamma > 1.) { // ?
			escf = 0.5*expint(2, d1);
		}
		else {	
			lvg_fc.td = d1;
			escf = 0.5*gamma *qromb<lvg_func_g>(lvg_fc, mu_c, 1., G_FUNC_PRECISION);
		}

		if (mu_c*mu_c*gamma > 1. && d2*d2*gamma > 1.) {
			escf += 0.5*expint(2, d2);
		}
		else {
			lvg_fc.td = d2;
			escf += 0.5*gamma *qromb<lvg_func_g>(lvg_fc, mu_c, 1., G_FUNC_PRECISION);
		}
	}
	return escf;
}

lvg_method_data_spline::lvg_method_data_spline(const std::string& path, std::string name, int verbosity)
    : lvg_method_data(path, name, verbosity)
{
    int i, j;
    p_d1_arr = alloc_2d_array<double>(nb_d, nb_g);
    p_d2_arr = alloc_2d_array<double>(nb_d, nb_g);
    p_d12_arr = alloc_2d_array<double>(nb_d, nb_g);

    // The estimation of the derivatives; derivative on delta (first variable):
    for (j = 0; j < nb_g; j++)
    {
        for (i = 1; i < nb_d - 1; i++)
        {
            p_d1_arr[i][j] = 0.5 * ((p_arr[i + 1][j] - p_arr[i][j]) / (delta_arr[i + 1] - delta_arr[i]) +
                (p_arr[i][j] - p_arr[i - 1][j]) / (delta_arr[i] - delta_arr[i - 1]));
        }
        p_d1_arr[0][j] = (p_arr[1][j] - p_arr[0][j]) / (delta_arr[1] - delta_arr[0]);
        p_d1_arr[nb_d - 1][j] = (p_arr[nb_d - 1][j] - p_arr[nb_d - 2][j]) / (delta_arr[nb_d - 1] - delta_arr[nb_d - 2]);
    }

    // derivative on gamma (second variable):
    for (i = 0; i < nb_d; i++)
    {
        for (j = 1; j < nb_g - 1; j++)
        {
            p_d2_arr[i][j] = 0.5 * ((p_arr[i][j + 1] - p_arr[i][j]) / (gamma_arr[j + 1] - gamma_arr[j]) +
                (p_arr[i][j] - p_arr[i][j - 1]) / (gamma_arr[j] - gamma_arr[j - 1]));
        }
        p_d2_arr[i][0] = (p_arr[i][1] - p_arr[i][0]) / (gamma_arr[1] - gamma_arr[0]);
        p_d2_arr[i][nb_g - 1] = (p_arr[i][nb_g - 1] - p_arr[i][nb_g - 2]) / (gamma_arr[nb_g - 1] - gamma_arr[nb_g - 2]);
    }

    // cross-derivative:
    for (i = 0; i < nb_d; i++)
    {
        for (j = 1; j < nb_g - 1; j++)
        {
            p_d12_arr[i][j] = 0.25 * ((p_d1_arr[i][j + 1] - p_d1_arr[i][j]) / (gamma_arr[j + 1] - gamma_arr[j]) +
                (p_d1_arr[i][j] - p_d1_arr[i][j - 1]) / (gamma_arr[j] - gamma_arr[j - 1]));
        }
        p_d12_arr[i][0] = 0.5 * (p_d1_arr[i][1] - p_d1_arr[i][0]) / (gamma_arr[1] - gamma_arr[0]);
        p_d12_arr[i][nb_g - 1] = 0.5 * (p_d1_arr[i][nb_g - 1] - p_d1_arr[i][nb_g - 2]) / (gamma_arr[nb_g - 1] - gamma_arr[nb_g - 2]);
    }

    for (j = 0; j < nb_g; j++)
    {
        for (i = 1; i < nb_d - 1; i++)
        {
            p_d12_arr[i][j] += 0.25 * ((p_d2_arr[i + 1][j] - p_d2_arr[i][j]) / (delta_arr[i + 1] - delta_arr[i]) +
                (p_d2_arr[i][j] - p_d2_arr[i - 1][j]) / (delta_arr[i] - delta_arr[i - 1]));
        }
        p_d12_arr[0][j] += 0.5 * (p_d2_arr[1][j] - p_d2_arr[0][j]) / (delta_arr[1] - delta_arr[0]);
        p_d12_arr[nb_d - 1][j] += 0.5 * (p_d2_arr[nb_d - 1][j] - p_d2_arr[nb_d - 2][j]) / (delta_arr[nb_d - 1] - delta_arr[nb_d - 2]);
    }
}

lvg_method_data_spline::~lvg_method_data_spline()
{
    free_2d_array(p_d1_arr);
    free_2d_array(p_d2_arr);
    free_2d_array(p_d12_arr);
}

double lvg_method_data_spline::get_esc_func(double gamma, double delta) const
{
    int k, l;
    double escf(0.), eg1, eg2, y_arr[4], y1_arr[4], y2_arr[4], y12_arr[4];

    // if values in question are out of range of the array with a dimension dim, returned index < 0 or index >= dim-1
    locate_index(delta_arr, nb_d, delta, k);
    locate_index(gamma_arr, nb_g, gamma, l);

    // in order to avoid wrong extrapolation:
    if (k < 0) {
        k = 0;
        delta = delta_arr[0];
    }
    else if (k > nb_d - 2) {
        k = nb_d - 2;
        delta = delta_arr[nb_d - 1];
    }

    if (l < 0) {
        l = 0;
        gamma = gamma_arr[0];
    }
    else if (l > nb_g - 2) {
        l = nb_g - 2;
        gamma = gamma_arr[nb_g - 1];
    }

    y_arr[0] = p_arr[k][l];
    y_arr[1] = p_arr[k + 1][l];
    y_arr[2] = p_arr[k + 1][l + 1];
    y_arr[3] = p_arr[k][l + 1];

    y1_arr[0] = p_d1_arr[k][l];
    y1_arr[1] = p_d1_arr[k + 1][l];
    y1_arr[2] = p_d1_arr[k + 1][l + 1];
    y1_arr[3] = p_d1_arr[k][l + 1];

    y2_arr[0] = p_d2_arr[k][l];
    y2_arr[1] = p_d2_arr[k + 1][l];
    y2_arr[2] = p_d2_arr[k + 1][l + 1];
    y2_arr[3] = p_d2_arr[k][l + 1];

    y12_arr[0] = p_d12_arr[k][l];
    y12_arr[1] = p_d12_arr[k + 1][l];
    y12_arr[2] = p_d12_arr[k + 1][l + 1];
    y12_arr[3] = p_d12_arr[k][l + 1];

    bcuint(y_arr, y1_arr, y2_arr, y12_arr, delta_arr[k], delta_arr[k + 1], gamma_arr[l], gamma_arr[l + 1], delta, gamma, escf, eg1, eg2);

    if (escf > 1.)
        escf = 1.;
    else if (escf < 0.)
        escf = 0.;

    return escf;
}

lvg_line_overlap_data::lvg_line_overlap_data(const std::string& path, std::string name, int verbosity)
{
    char text_line[MAX_TEXT_LINE_WIDTH];
    int i, j, k, l;

    string file_name;
    ifstream input;

    file_name = path + name;
    input.open(file_name.c_str(), ios_base::in);

    if (!input.is_open()) {
        cout << "Error in lvg_method_functions.cpp: can't open " << file_name << endl;
        exit(1);
    }

    // comment lines
    input.getline(text_line, MAX_TEXT_LINE_WIDTH);
    input.getline(text_line, MAX_TEXT_LINE_WIDTH);
    input.getline(text_line, MAX_TEXT_LINE_WIDTH);
    input >> nb_d >> nb_dx >> nb_gr >> nb_g;

    delta_arr = new double[nb_d];
    gamma_arr = new double[nb_g];
    gratio_arr = new double[nb_gr];
    dx_arr = new double[nb_dx];

    p_arr = alloc_2d_array<double>(nb_dx*nb_d, nb_g*nb_gr);
    for (i = 0; i < nb_g; i++) {
        input >> gamma_arr[i];
    }

    for (l = 0; l < nb_d; l++) {
        for (k = 0; k < nb_dx; k++) {
            input >> delta_arr[l] >> dx_arr[k];
            delta_arr[l] = log10(delta_arr[l]); // logarithmic scale on delta

            for (i = 0; i < nb_gr; i++) {
                input >> j >> gratio_arr[i];
                for (j = 0; j < nb_g; j++) {
                    input >> p_arr[l * nb_dx + k][i * nb_g + j];
                    // Note!!! the current version of lvg data for second line are normalized incorrectly,
                    if (name == "lvg/line_overlap_func_p2_test.txt") {
                        p_arr[l * nb_dx + k][i * nb_g + j] /= gratio_arr[i];
                    }
                }
            }
        }
    }
    input.close();

    if (verbosity)
        cout << "The data on LVG method (line overlaps) have been initialized." << endl;
}

lvg_line_overlap_data::~lvg_line_overlap_data()
{
    delete[] delta_arr;
    delete[] gamma_arr;
    delete[] gratio_arr;
    delete[] dx_arr;
    free_2d_array(p_arr);
}

double lvg_line_overlap_data::get_esc_func(double gamma, double delta, double gamma_ratio, double delta_x) const
{
    int l, k, n, m;
    double u, t, p, y, escf(0.);

    delta = log10(delta);
    locate_index(delta_arr, nb_d, delta, m);
    locate_index(gamma_arr, nb_g, gamma, l);    
    locate_index(gratio_arr, nb_gr, gamma_ratio, k);
    locate_index(dx_arr, nb_dx, delta_x, n);

    y = u = t = p = 0.;

    if (m < 0)
        m = 0;
    else if (m > nb_d - 2) {
        m = nb_d - 2;
        y = 1.;
    }
    else y = (delta - delta_arr[m]) / (delta_arr[m + 1] - delta_arr[m]);

    if (n < 0)
        n = 0;
    else if (n > nb_dx - 2) {
        p = 1.;
        n = nb_dx - 2;
    }
    else p = (delta_x - dx_arr[n]) / (dx_arr[n + 1] - dx_arr[n]);

    if (l < 0) 
        l = 0;
    else if (l > nb_g - 2) {
        l = nb_g - 2;
        u = 1.;
    }
    else u = (gamma - gamma_arr[l]) / (gamma_arr[l + 1] - gamma_arr[l]);

    if (k < 0) 
        k = 0;
    else if (k > nb_gr - 2) {
        t = 1.;
        k = nb_gr - 2;
    }
    else t = (gamma_ratio - gratio_arr[k]) / (gratio_arr[k + 1] - gratio_arr[k]);

    escf = p_arr[m * nb_dx + n][k * nb_g + l] * (1. - u) * (1. - t) * (1. - p) * (1. - y)
        + p_arr[m * nb_dx + n][k * nb_g + l + 1] * u *(1. - t) * (1. - p) * (1. - y)
        + p_arr[m * nb_dx + n][(k + 1) * nb_g + l] * (1. - u) * t* (1. - p)* (1. - y)
        + p_arr[m * nb_dx + n][(k + 1) * nb_g + l + 1] * u * t * (1. - p) * (1. - y)
        + p_arr[m * nb_dx + n + 1][k * nb_g + l] * (1. - u) * (1. - t)* p * (1. - y)
        + p_arr[m * nb_dx + n + 1][k * nb_g + l + 1] * u * (1. - t) * p * (1. - y)
        + p_arr[m * nb_dx + n + 1][(k + 1) * nb_g + l] * (1. - u) * t * p * (1. - y)
        + p_arr[m * nb_dx + n + 1][(k + 1) * nb_g + l + 1] * u * t * p * (1. - y)
        + p_arr[(m + 1) * nb_dx + n][k * nb_g + l] * (1. - u) * (1. - t) * (1. - p) * y
        + p_arr[(m + 1) * nb_dx + n][k * nb_g + l + 1] * u * (1. - t) * (1. - p) * y
        + p_arr[(m + 1) * nb_dx + n][(k + 1) * nb_g + l] * (1. - u) *t * (1. - p) * y
        + p_arr[(m + 1) * nb_dx + n][(k + 1) * nb_g + l + 1] * u * t * (1. - p)* y
        + p_arr[(m + 1) * nb_dx + n + 1][k * nb_g + l] * (1. - u) * (1. - t) * p * y
        + p_arr[(m + 1) * nb_dx + n + 1][k * nb_g + l + 1] * u * (1. - t) * p * y
        + p_arr[(m + 1) * nb_dx + n + 1][(k + 1) * nb_g + l] * (1. - u) *t *p * y 
        + p_arr[(m + 1) * nb_dx + n + 1][(k + 1) * nb_g + l + 1] * u* t* p* y;

    if (escf > 1.)
        escf = 1.;
    else if (escf < 0.)
        escf = 0.;

    return escf;
}
