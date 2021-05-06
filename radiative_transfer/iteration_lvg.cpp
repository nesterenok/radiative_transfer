
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <omp.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cfloat>
#include <cstring>
#include <stdlib.h>

#include "utils.h"
#include "linear_algebra.h"
#include "special_functions.h"
#include "integration.h"
#include "interpolation.h"
#include "constants.h"
#include "iteration_lvg.h"

#define SOURCE_NAME "iteration_lvg.cpp"
#define INV_TRANS_FACTOR -0.1
#define MIN_LINE_OPACITY 1.e-99  // arbitrary very small value;
using namespace std;

// The molecule data are not initialized in the constructor;
iteration_scheme_lvg::iteration_scheme_lvg(const std::string& data_path, int v)
	: diagram(0), einst_coeff(0), coll_trans(0), dust(0), 
    verbosity(v), nb_mol_lev(0), temp_n(0.), temp_el(0.), vel_grad(0.), vel_width(0.), mol_conc(0.),
	coll_partn_conc(0), indices(0), matrix(0)
{
	// LVG data:
	loss_func_line_phot = new lvg_method_data(data_path, "lvg/lvg_loss_func.txt", verbosity);
}

iteration_scheme_lvg::~iteration_scheme_lvg() {
    delete loss_func_line_phot;
    delete[] indices;
    delete[] coll_partn_conc;
	if (matrix != 0) 
        free_2d_array<double>(matrix);
}

void iteration_scheme_lvg::init_molecule_data(const energy_diagram *e_d, const einstein_coeff *e_c, const collisional_transitions *c_t)
{
	diagram = e_d;
	einst_coeff = e_c; 
	coll_trans = c_t;
	nb_mol_lev = diagram->nb_lev;

	// initialization of transition rate array:
	if (matrix != 0) 
		free_2d_array<double>(matrix);
	matrix = alloc_2d_array<double>(nb_mol_lev, nb_mol_lev);
}

// The molecular data must be initialized before the call of this function;
void iteration_scheme_lvg::set_parameters(double tn, double te, double el_conc, double h_conc, double ph2_conc, double oh2_conc, 
	double he_conc, double mc, double vel_turb)
{
	temp_n = tn;
	temp_el = te;
	mol_conc = mc;
	vel_width = pow(2.*BOLTZMANN_CONSTANT *temp_n /diagram->mol.mass + vel_turb *vel_turb, 0.5);		

	coll_trans->set_gas_param(temp_n, temp_el, he_conc, ph2_conc, oh2_conc, h_conc, el_conc, coll_partn_conc, indices);
}

void iteration_scheme_lvg::set_dust_parameters(vector<double>& conc, vector<double>& temp)
{
	dgrain_conc.clear();
	dgrain_temp.clear();

    for (int i = 0; i < (int)conc.size(); i++) {
        dgrain_conc.push_back(conc[i]);
    }
    for (int i = 0; i < (int)temp.size(); i++) {
		dgrain_temp.push_back(temp[i]);
	}
	if (dust->nb_of_comp != (int)dgrain_conc.size() || dust->nb_of_comp != (int)dgrain_temp.size()) {
		cout << "Error in " << SOURCE_NAME << ": there is an inconsistency in the interstellar dust data" << endl;
		exit(1);
	}
}

void iteration_scheme_lvg::calc_new_pop(double *old_pop, double *new_pop, double &eq_error)
{
	int i;
	double *f_vector, *b_vector;

	f_vector = new double [nb_mol_lev];
	b_vector = new double [nb_mol_lev];

	(*this)(nb_mol_lev, old_pop, f_vector);

	memset(b_vector, 0, nb_mol_lev*sizeof(double));
	b_vector[0] = 1.;

	lu_matrix_solve(matrix, b_vector, nb_mol_lev);
	memcpy(new_pop, b_vector, nb_mol_lev*sizeof(double));
	
	eq_error = 0.;
	for (i = 0; i < nb_mol_lev; i++) {
		if (eq_error < fabs(f_vector[i])) 
            eq_error = fabs(f_vector[i]);
	}
	delete [] f_vector;
	delete [] b_vector;
}

void iteration_scheme_lvg::operator() (int dim, double *level_pop, double *df)
{
	int i, j;
	double a, b, y, down_rate, up_rate, intensity;
	
	memset(*matrix, 0, nb_mol_lev*nb_mol_lev*sizeof(double));
	for (i = 1; i < nb_mol_lev; i++) {
		for (j = 0; j < i; j++)
		{
			coll_trans->get_rate_neutrals(diagram->lev_array[i], diagram->lev_array[j], down_rate, up_rate, temp_n, coll_partn_conc, indices);
			coll_trans->get_rate_electrons(diagram->lev_array[i], diagram->lev_array[j], a, b, temp_el, coll_partn_conc, indices);
				
			down_rate += a;
			up_rate += b;

			matrix[i][i] -= down_rate;
			matrix[j][i] += down_rate;

			matrix[i][j] += up_rate;
			matrix[j][j] -= up_rate;

			//if ((level_pop[j] *einst_coeff->arr[j][i] - level_pop[i] *einst_coeff->arr[i][j]) != 0.) 
			if (einst_coeff->arr[i][j] > 0.)
			{
				intensity_calc(i, j, level_pop, intensity);
				y = einst_coeff->arr[i][j] * (1. + intensity);

				matrix[i][i] -= y;
				matrix[j][i] += y;

				y = einst_coeff->arr[j][i] * intensity;
				matrix[i][j] += y;
				matrix[j][j] -= y;
			}
		}
	}

	for (j = 0; j < nb_mol_lev; j++) {
		matrix[0][j] = 1.;
	}
	
	memset(df, 0, nb_mol_lev*sizeof(double));
	df[0] = 1.;

	for (i = 0; i < nb_mol_lev; i++) {
		for (j = 0; j < nb_mol_lev; j++) {
			df[i] -= matrix[i][j] *level_pop[j];
		}
	}
}

void iteration_scheme_lvg::intensity_calc(int up, int low, double *level_pop, double &intensity) const
{
	double c, energy, dust_opacity, line_emiss, line_opacity, gamma, delta, ep1;
	
	energy = diagram->lev_array[up].energy - diagram->lev_array[low].energy;
	c = mol_conc/(EIGHT_PI * vel_width *energy*energy*energy);
	
    // line absorption and emission coefficients without 1/sqrt(pi),
	line_emiss = c*einst_coeff->arr[up][low] *level_pop[up];
	line_opacity = c*einst_coeff->arr[low][up] *level_pop[low] - line_emiss + MIN_LINE_OPACITY;
	
	if (line_opacity < 0.) 
		line_opacity *= INV_TRANS_FACTOR;

	dust_opacity = dust->absorption(energy, dgrain_conc);
	
    // there is no distinction between the cases with different dv/dz signs,
	gamma = fabs(vel_grad)/(vel_width*line_opacity);
	delta = fabs(vel_grad)/(vel_width*dust_opacity);

	ep1 = loss_func_line_phot->get_esc_func(gamma, delta);
	intensity = line_emiss/line_opacity *ep1; // without dust emission
};

void iteration_scheme_lvg::calc_line_stat(const string& fname, double* lev_pop)
{
    int i, j;
    double c, line_opacity, dust_opacity, energy, g, d;
    
    line_parameters lp;
    vector<line_parameters> lp_vector;
    ofstream output;

    for (i = 1; i < nb_mol_lev; i++) {
        for (j = 0; j < i; j++) {
            // the lines without population inversion are considered;
            if ((lev_pop[j] * einst_coeff->arr[j][i] - lev_pop[i] * einst_coeff->arr[i][j]) > 0.) {
                energy = diagram->lev_array[i].energy - diagram->lev_array[j].energy;
                c = mol_conc / (EIGHT_PI * vel_width * energy * energy * energy);
                
                dust_opacity = dust->absorption(energy, dgrain_conc);
                line_opacity = c * (einst_coeff->arr[j][i] * lev_pop[j] - einst_coeff->arr[i][j] * lev_pop[i]) + MIN_LINE_OPACITY;

                g = fabs(vel_grad) / (vel_width * line_opacity); // absolute value
                d = fabs(vel_grad) / (vel_width * dust_opacity);

                lp.nbl = j;
                lp.nbu = i;
                lp.g = g;
                lp.d = d;
                lp.en = energy;
                lp_vector.push_back(lp);
            }
        }
    }
 
    output.open(fname.c_str(), std::ios_base::out);
    if (!output.is_open())
        cout << "Error in " << SOURCE_NAME << ": can't open file to write line statistics data;" << endl;
    else {
        output << scientific;
        output.precision(4);

        output << "! parameters: upper, lower levels, energy(cm-1), abs(d), abs(g)" << endl;
        for (i = 0; i < (int)lp_vector.size(); i++) {
            output << left << setw(5) << lp_vector[i].nbu << setw(5) << lp_vector[i].nbl << setw(15) << lp_vector[i].en
                << setw(15) << lp_vector[i].d << setw(15) << lp_vector[i].g << endl;
        }
        output.close();
    }
}


hfs_lines::hfs_lines() : nb(0)
{;}

void hfs_lines::clear() {
    nb = 0;
    upl.clear();
    lowl.clear();
    en.clear();
}

void hfs_lines::add_line(int u, int l, double e) {
    upl.push_back(u);
    lowl.push_back(l);
    en.push_back(e);
    nb++;
}

void hfs_lines::exchange_lines(int i, int j) {
    swap(upl[i], upl[j]);
    swap(lowl[i], lowl[j]);
    swap(en[i], en[j]);
}

void hfs_lines::sort()
{
    if (nb <= 2)
        return;

    int i, j;
    double e;
    for (i = 0; i < nb; i++) {
        for (j = i + 1; j < nb; j++) {
            if (en[j] < en[i]) {
                exchange_lines(i, j);
            }
        }
    }
    // finding the lines most closest in frequency,
    j = 0;
    e = en[1] - en[0];
    for (i = 1; i < nb - 1; i++) {
        if (en[i + 1] - en[i] < e) {
            j = i;
        }
    }
    e = en[j];

    for (i = 0; i < nb; i++) {
        for (j = i + 1; j < nb; j++) {
            if (fabs(en[j] - e) < fabs(en[i] - e)) {
                exchange_lines(i, j);
            }
        }
    }
}

void hfs_lines::split(hfs_lines & hfs_l) 
{
    hfs_l.clear();
    hfs_l.add_line(upl[0], lowl[0], en[0]);
    hfs_l.add_line(upl[1], lowl[1], en[1]);

    upl.erase(upl.begin(), upl.begin() + 2);
    lowl.erase(lowl.begin(), lowl.begin() + 2);
    en.erase(en.begin(), en.begin() + 2);
    nb -= 2;
}


iteration_scheme_line_overlap::iteration_scheme_line_overlap(const std::string& data_path, int verbosity)
    :iteration_scheme_lvg(data_path, verbosity) {
    max_dx = 4; // > 0., may differ from the max value given in data files
    line_overlap1 = new lvg_line_overlap_data(data_path, "lvg/line_overlap_func_p1.txt", verbosity);
    line_overlap2 = new lvg_line_overlap_data(data_path, "lvg/line_overlap_func_p2.txt", verbosity);
}

iteration_scheme_line_overlap::~iteration_scheme_line_overlap() {
    delete line_overlap1;
    delete line_overlap2;
}

void iteration_scheme_line_overlap::init_molecule_data(const energy_diagram* e_d, const einstein_coeff* e_c, const collisional_transitions* c_t)
{
    int i, j, m, l;
    hfs_lines hfs_l, hfs_l2;
    
    iteration_scheme_lvg::init_molecule_data(e_d, e_c, c_t);
    line_list.clear();

    for (i = 2; i < nb_mol_lev; i += 2) {
        for (j = 0; j < i; j += 2)
        {
            hfs_l.clear();
            for (m = 0; m < 2; m++) {
                for (l = 0; l < 2; l++) {
                    if (einst_coeff->arr[i + m][j + l] > 1.e-99) {
                        hfs_l.add_line(i + m, j + l, diagram->lev_array[i + m].energy - diagram->lev_array[j + l].energy);
                    }
                }
            }          
            hfs_l.sort();

            if (hfs_l.nb >= 3) { // nb may be 0,1,2,3,4
                hfs_l.split(hfs_l2);
                line_list.push_back(hfs_l2);
            }
            if (hfs_l.nb > 0)
                line_list.push_back(hfs_l);
        }
    }
}

void iteration_scheme_line_overlap::operator()(int dim, double* level_pop, double* df)
{
    int i, j, m, l, mm, ll;
    double a, b, y, down_rate, up_rate, intens1, intens2;

    memset(*matrix, 0, nb_mol_lev * nb_mol_lev * sizeof(double));
    for (i = 1; i < nb_mol_lev; i++) {
        for (j = 0; j < i; j++)
        {
            coll_trans->get_rate_neutrals(diagram->lev_array[i], diagram->lev_array[j], down_rate, up_rate, temp_n, coll_partn_conc, indices);
            coll_trans->get_rate_electrons(diagram->lev_array[i], diagram->lev_array[j], a, b, temp_el, coll_partn_conc, indices);

            down_rate += a;
            up_rate += b;

            matrix[i][i] -= down_rate;
            matrix[j][i] += down_rate;

            matrix[i][j] += up_rate;
            matrix[j][j] -= up_rate;
        }
    }

    for (i = 0; i < (int) line_list.size(); i++) {
        if (line_list[i].nb == 1) 
        {
            m = line_list[i].upl[0];
            l = line_list[i].lowl[0];

            iteration_scheme_lvg::intensity_calc(m, l, level_pop, intens1);
            y = einst_coeff->arr[m][l] * (1. + intens1);

            matrix[m][m] -= y;
            matrix[l][m] += y;

            y = einst_coeff->arr[l][m] * intens1;
            matrix[m][l] += y;
            matrix[l][l] -= y;
        }
        else if (line_list[i].nb == 2) {
            m = line_list[i].upl[0];
            l = line_list[i].lowl[0];

            mm = line_list[i].upl[1];
            ll = line_list[i].lowl[1];
            
            intensity_calc(m, l, mm, ll, level_pop, intens1, intens2);
            
            y = einst_coeff->arr[m][l] * (1. + intens1);
            matrix[m][m] -= y;
            matrix[l][m] += y;

            y = einst_coeff->arr[l][m] * intens1;
            matrix[m][l] += y;
            matrix[l][l] -= y;

            y = einst_coeff->arr[mm][ll] * (1. + intens2);
            matrix[mm][mm] -= y;
            matrix[ll][mm] += y;

            y = einst_coeff->arr[ll][mm] * intens2;
            matrix[mm][ll] += y;
            matrix[ll][ll] -= y;
        }
    }

    for (j = 0; j < nb_mol_lev; j++) {
        matrix[0][j] = 1.;
    }

    memset(df, 0, nb_mol_lev * sizeof(double));
    df[0] = 1.;

    for (i = 0; i < nb_mol_lev; i++) {
        for (j = 0; j < nb_mol_lev; j++) {
            df[i] -= matrix[i][j] * level_pop[j];
        }
    }
}

void iteration_scheme_line_overlap::intensity_calc(int u1, int l1, int u2, int l2, double* level_pop, double& intens1, double& intens2) const
{
    double c, energy, line_emiss1, line_opacity1, line_emiss2, line_opacity2, gamma1, gamma2, gratio, dx, ep1, ep2, ep01, ep02,
        dust_opacity, delta;

    energy = diagram->lev_array[u1].energy - diagram->lev_array[l1].energy;
    c = mol_conc / (EIGHT_PI * vel_width * energy * energy * energy);

    line_emiss1 = c * einst_coeff->arr[u1][l1] * level_pop[u1];
    line_opacity1 = c * (einst_coeff->arr[l1][u1] * level_pop[l1] - einst_coeff->arr[u1][l1] * level_pop[u1]) + MIN_LINE_OPACITY;

    line_emiss2 = c * einst_coeff->arr[u2][l2] * level_pop[u2];
    line_opacity2 = c * (einst_coeff->arr[l2][u2] * level_pop[l2] - einst_coeff->arr[u2][l2] * level_pop[u2]) + MIN_LINE_OPACITY;

    if (line_opacity1 < 0.)
        line_opacity1 *= INV_TRANS_FACTOR;

    if (line_opacity2 < 0.)
        line_opacity2 *= INV_TRANS_FACTOR;
 
    gamma1 = fabs(vel_grad) / (vel_width * line_opacity1);
    gamma2 = fabs(vel_grad) / (vel_width * line_opacity2);

    dust_opacity = dust->absorption(energy, dgrain_conc);
    delta = fabs(vel_grad) / (vel_width * dust_opacity); // is common parameter for both lines,

    dx = (diagram->lev_array[u1].energy - diagram->lev_array[l1].energy - diagram->lev_array[u2].energy + diagram->lev_array[l2].energy)
       * SPEED_OF_LIGHT / (energy * vel_width);
    
    if (vel_grad < 0.)
        dx *= -1.;

    ep1 = ep2 = ep01 = ep02 = 0.;
    if (fabs(dx) < max_dx) {
        gratio = gamma2 / gamma1;
        ep1 = line_overlap1->get_esc_func(gamma1, delta, gratio, dx);

        gratio = gamma1 / gamma2;
        ep2 = line_overlap1->get_esc_func(gamma2, delta, gratio, -dx);
    }

    if (fabs(dx) > max_dx - 0.5) {
        ep01 = loss_func_line_phot->get_esc_func(gamma1, delta);
        ep02 = loss_func_line_phot->get_esc_func(gamma2, delta);   
    }
 
    if (fabs(dx) > max_dx) {
         ep1 = ep01;
         ep2 = ep02;
    }
    else if (fabs(dx) > max_dx - 0.5) {
        c = 2.*(max_dx - fabs(dx));
        ep1 = ep01 * (1. - c) + ep1 * c;
        ep2 = ep02 * (1. - c) + ep2 * c;
    }
    intens1 = line_emiss1 / line_opacity1 * ep1;
    intens2 = line_emiss2 / line_opacity2 * ep2;

    if (fabs(dx) < max_dx) {
        gratio = gamma2 / gamma1;
        ep1 = line_overlap2->get_esc_func(gamma1, delta, gratio, dx);

        gratio = gamma1 / gamma2;
        ep2 = line_overlap2->get_esc_func(gamma2, delta, gratio, -dx);

        if (fabs(dx) > max_dx - 0.5) {
            c = 2.*(max_dx - fabs(dx));
            ep1 *= c;
            ep2 *= c;
        }
        intens1 += line_emiss2 / line_opacity2 * ep1;
        intens2 += line_emiss1 / line_opacity1 * ep2;
    }
}

void iteration_scheme_line_overlap::intensity_calc(int upl, int lowl, double* level_pop, double& intensity) const
{
    int i, k;
    double a;

    intensity = 0.;
    for (i = 0; i < (int) line_list.size(); i++) {
        for (k = 0; k < line_list[i].nb; k++) 
        {
            if (line_list[i].lowl[k] == lowl && line_list[i].upl[k] == upl) {
                if (line_list[i].nb == 1) {
                    iteration_scheme_lvg::intensity_calc(upl, lowl, level_pop, intensity);
                }
                else { // probably, the order is not important 
                    if (k == 0) 
                        intensity_calc(upl, lowl, line_list[i].upl[1], line_list[i].lowl[1], level_pop, intensity, a);
                    else if (k == 1)
                        intensity_calc(line_list[i].upl[0], line_list[i].lowl[0], upl, lowl, level_pop, a, intensity);
                }
                return;
            }
        }
    }
}

void get_nb_overlap_lines(const energy_diagram* diagram, const einstein_coeff* einst_coeff, double vel_width, 
    int& nb_double_overlap, int& nb_triple_overlap)
{
    const double dx_lim = 4.;
    int i, j, m, l, nb_mol_lev;
    double en, dx;   
    hfs_lines hfs_l, hfs_l2;
 
    nb_mol_lev = diagram->nb_lev;
    nb_double_overlap = nb_triple_overlap = 0;

    for (i = 2; i < nb_mol_lev; i += 2) {
        for (j = 0; j < i; j += 2)
        {
            hfs_l.clear();
            for (m = 0; m < 2; m++) {
                for (l = 0; l < 2; l++) {
                    if (einst_coeff->arr[i + m][j + l] > 1.e-99) {
                        hfs_l.add_line(i + m, j + l, diagram->lev_array[i + m].energy - diagram->lev_array[j + l].energy);
                    }
                }
            }
            hfs_l.sort();

            // it is assumed that line frequency >> frequency difference between lines,
            en = diagram->lev_array[i].energy - diagram->lev_array[j].energy;
            if (hfs_l.nb > 1) {
                dx = (diagram->lev_array[hfs_l.upl[0]].energy - diagram->lev_array[hfs_l.lowl[0]].energy
                    - diagram->lev_array[hfs_l.upl[1]].energy + diagram->lev_array[hfs_l.lowl[1]].energy)
                    * SPEED_OF_LIGHT / (en * vel_width);

                if (fabs(dx) < dx_lim)
                    nb_double_overlap++;
            }
            if (hfs_l.nb > 2) {
                dx = (diagram->lev_array[hfs_l.upl[1]].energy - diagram->lev_array[hfs_l.lowl[1]].energy
                    - diagram->lev_array[hfs_l.upl[2]].energy + diagram->lev_array[hfs_l.lowl[2]].energy)
                    * SPEED_OF_LIGHT / (en * vel_width);

                if (fabs(dx) < dx_lim)
                    nb_triple_overlap++;
            }

            if (hfs_l.nb == 4) {
                hfs_l.split(hfs_l2);
                
                dx = (diagram->lev_array[hfs_l.upl[0]].energy - diagram->lev_array[hfs_l.lowl[0]].energy
                    - diagram->lev_array[hfs_l.upl[1]].energy + diagram->lev_array[hfs_l.lowl[1]].energy)
                    * SPEED_OF_LIGHT / (en * vel_width);

                if (fabs(dx) < dx_lim)
                    nb_double_overlap++;
            }
        }
    }
}

/*
    //void calc_line_stat(const std::string& path, const std::string& str_id, double*);

struct line_param {
    int v_l, v_u, nb_l, nb_u;
    double g, d, ep0, ep1, ep2, sc, intens;
};


void iteration_scheme_lvg::calc_line_stat(const string &path, const std::string &str_id, double *level_pop)
{
    int i, j, nb, nb2, nb_points_bin;
    double line_opacity, dust_opacity, energy, min, max, p, dust_emiss;
    double *x_arr, *y_arr1, *y_arr2, *y_arr3, *g_arr;

    string fname;
    ofstream output;

    line_param lp;
    vector<line_param> line_param_vect;

    lvg_func_slab calc;

    for (i = 1; i < nb_mol_lev; i++) {
        for (j = 0; j < i; j++)
        {
            // The lines without population inversion are considered;
            if ((level_pop[j] *einst_coeff->arr[j][i] - level_pop[i] *einst_coeff->arr[i][j]) > 0.)
            {
                lp.v_l = diagram->lev_array[j].v;
                lp.v_u = diagram->lev_array[i].v;

                lp.nb_l = diagram->lev_array[j].nb;
                lp.nb_u = diagram->lev_array[i].nb;

                energy = diagram->lev_array[i].energy - diagram->lev_array[j].energy;

                dust_emiss = dust->emissivity(energy, dust_temp) *dust_density;
                dust_opacity = dust->opacity(energy) *dust_density;
                line_opacity = abs_data *(einst_coeff->arr[j][i] *level_pop[j] - einst_coeff->arr[i][j] *level_pop[i]) /(energy*energy*energy);

                lp.g = vel_grad/(vel_width *line_opacity);
                lp.d = vel_grad/(vel_width *dust_opacity);
                lp.sc = dust_emiss/dust_opacity;

                calc.g = lp.g;
                lp.ep0 = lp.g*qromb(calc, 1.e-6, 1., 1.e-6);

                if (lp.ep0 > 1.) lp.ep0 = 1.;
                else if (lp.ep0 < 0.) lp.ep0 = 0.;

                intensity_calc(i, j, level_pop);

                lp.ep1 = escp1;
                lp.ep2 = escp2;
                lp.intens = intensity;

                if (lp.ep1 > 1.-lp.ep0) {
                    lp.ep1 = 1.-lp.ep0;
                }
                line_param_vect.push_back(lp);
            }
        }
    }
    nb = init_log_grid(x_arr, min =1.e-8, max =1., nb_points_bin =5);
    y_arr1 = new double [nb];
    y_arr2 = new double [nb];

    memset(y_arr1, 0, nb*sizeof(double));
    memset(y_arr2, 0, nb*sizeof(double));

    nb2 = init_log_grid(g_arr, min =1.e-8, max =1.e+16, nb_points_bin =5);
    y_arr3 = new double [nb2];
    memset(y_arr3, 0, nb2*sizeof(double));

    for (i = 0; i < (int) line_param_vect.size(); i++)
    {
        p = line_param_vect[i].ep1/(1. - line_param_vect[i].ep0 + DBL_EPSILON);
        locate_index(x_arr, nb, p, j);
        if (j >= 0 && j < nb-1) y_arr1[j] += 1.;

        p = line_param_vect[i].sc*(1. - line_param_vect[i].ep1 - line_param_vect[i].ep2) /(line_param_vect[i].intens + DBL_EPSILON);
        locate_index(x_arr, nb, p, j);
        if (j >= 0 && j < nb-1) y_arr2[j] += 1.;

        locate_index(g_arr, nb2, line_param_vect[i].g, j);
        if (j >= 0 && j < nb2-1) y_arr3[j] += 1.;
    }

    fname = path + diagram->mol.name;
    fname += "_line_stat";
    fname += str_id;
    fname += "_1.txt";
    output.open(fname.c_str(), std::ios_base::out);
    output << scientific;
    output.precision(4);

    output << "Nb of lines: " << line_param_vect.size() << endl;
    output << "The distribution of the parameter g:" << endl;
    for (i = 0; i < nb2-1; i++)
    {
        y_arr3[i] /= (double) line_param_vect.size();
        output << left << setw(15) << g_arr[i] << setw(15) << y_arr3[i] << endl
            << setw(15) << g_arr[i+1] << setw(15) << y_arr3[i] << endl;
    }

    output << endl << "The impact of dust to the escape probability in the line source function term;" << endl;
    for (i = 0; i < nb-1; i++)
    {
        y_arr1[i] /= (double) line_param_vect.size();
        output << left << setw(15) << x_arr[i] << setw(15) << y_arr1[i] << endl
            << setw(15) << x_arr[i+1] << setw(15) << y_arr1[i] << endl;
    }

    output << endl << "The impact of dust term to the intensity;" << endl;
    for (i = 0; i < nb-1; i++)
    {
        y_arr2[i] /= (double) line_param_vect.size();
        output << left << setw(15) << x_arr[i] << setw(15) << y_arr2[i] << endl
            << setw(15) << x_arr[i+1] << setw(15) << y_arr2[i] << endl;
    }
    output.close();

    fname = path + diagram->mol.name;
    fname += "_line_stat";
    fname += str_id;
    fname += "_2.txt";
    output.open(fname.c_str(), std::ios_base::out);
    output << scientific;
    output.precision(4);

    output << "The molecule lines without change of the vibrational state;" << endl;
    output << "d, g" << endl;
    for (i = 0; i < (int) line_param_vect.size(); i += 10)
    {
        if (line_param_vect[i].v_l == line_param_vect[i].v_u) {
            output << left << setw(5) << line_param_vect[i].nb_l << setw(5) << line_param_vect[i].nb_u
                << setw(15) << line_param_vect[i].d << setw(15) << line_param_vect[i].g << endl;
        }
    }
    output.close();

    fname = path + diagram->mol.name;
    fname += "_line_stat";
    fname += str_id;
    fname += "_3.txt";
    output.open(fname.c_str(), std::ios_base::out);
    output << scientific;
    output.precision(4);

    output << "The molecule lines with change of the vibrational state;" << endl;
    output << "d, g" << endl;
    for (i = 0; i < (int) line_param_vect.size(); i += 10)
    {
        if (line_param_vect[i].v_l != line_param_vect[i].v_u) {
            output << left << setw(5) << line_param_vect[i].nb_l << setw(5) << line_param_vect[i].nb_u
                << setw(15) << line_param_vect[i].d << setw(15) << line_param_vect[i].g << endl;
        }
    }
    output.close();

    delete [] x_arr;
    delete [] y_arr1;
    delete [] y_arr2;
    delete [] y_arr3;
    delete [] g_arr;
}
*/
