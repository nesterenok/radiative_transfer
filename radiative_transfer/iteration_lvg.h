#pragma once
#include <math.h>
#include <vector>
#include "constants.h"
#include "spectroscopy.h"
#include "coll_rates.h"
#include "lvg_method_functions.h"
#include "dust_model.h"

class line_parameters {
public:
    int nbl, nbu;
    double g, d, en;
};

// it is necessary to add the capability of reducing the number of levels during calculations,
class iteration_scheme_lvg
{
protected:
	int		nb_mol_lev, verbosity;
	double	temp_n, temp_el, mol_conc, vel_width, vel_grad;
    
    int		*indices;
	double	*coll_partn_conc, **matrix;
	std::vector<double> dgrain_conc, dgrain_temp;
	
	const energy_diagram	*diagram;
	const einstein_coeff	*einst_coeff;
	const collisional_transitions *coll_trans;
    
    const dust_model		*dust;	
    const lvg_method_data   *loss_func_line_phot;
	std::vector<const radiation_field*> ext_rad_field;

public:
	int get_vector_dim() { return nb_mol_lev; } // is used in iteration_control class
	
	// these function has to be called if molecule is changed:
	virtual void init_molecule_data(const energy_diagram *, const einstein_coeff *, const collisional_transitions *);

	void set_parameters(double temp_n, double temp_e, double el_conc, double h_conc, double ph2_conc, double oh2_conc, double he_conc, 
		double mol_conc, double vel_turb);
	void set_dust_parameters(std::vector<double> & conc, std::vector<double> & temp);
	void set_vel_grad(double vg) { vel_grad = vg; }
	
	void add_radiation_field(const radiation_field *erf) { ext_rad_field.push_back(erf); }
	void set_dust_model(const dust_model *d) { dust = d; }

    // intensity is dimensionless, = lambda^2/(2hv)*I
	virtual void intensity_calc(int upl, int lowl, double *lev_pop, double &intensity) const;
    
    // eq_error is the maximal error of equations for level populations, max df[], based on old populations
	void calc_new_pop(double* old_pop, double* new_pop, double& eq_error);
	
	// this function is necessary for using this class in Newton-Raphson algorithm;
	// dim is the dimension of the array, df[] - right-hand side residuals f_given[] - f_calc[]
	virtual void operator() (int dim, double *level_pop, double *df);
	
    // calculation and saving the spectral line parameters,
    void calc_line_stat(const std::string& fname, double *lev_pop);

    // path to the input data directory,
	iteration_scheme_lvg(const std::string & data_path, int verbosity = 1);
	virtual ~iteration_scheme_lvg();
};

// Hyper-fine splitted lines,
// lines that have levels with equal quantum numbers except the total (with HF splitting) angular momentum,
class hfs_lines 
{
public:
    int nb; // nb of lines in the set
    std::vector<int> upl, lowl;
    std::vector<double> en; // line energy in cm-1

    void add_line(int upl, int lowl, double energy);
    void exchange_lines(int i, int j);

    void sort();
    void split(hfs_lines & hfs_l); // takes away first two and creates from them new object
    void clear();
    hfs_lines();
};

// Only is suitable  for OH! The finding of overlapping lines uses explicitly the energy level structure;
// the maximal number of overlapped lines is equal to 2,
class iteration_scheme_line_overlap : public iteration_scheme_lvg
{
protected:
    double max_dx; // may be different tha the max value in file data;
    // lines are grouped in the sets of HF splitted lines, the list contains all lines, the number of lines in one group - 1 or 2; 
    std::vector<hfs_lines> line_list; 
    const lvg_line_overlap_data *line_overlap1, *line_overlap2;

public:
    void init_molecule_data(const energy_diagram*, const einstein_coeff*, const collisional_transitions*);
    
    void operator() (int dim, double* level_pop, double* df);
    void intensity_calc(int u1, int l1, int u2, int l2, double* level_pop, double& intens1, double& intens2) const;

    // the function is needed for limiting luminosity calculations,
    void intensity_calc(int upl, int lowl, double *lev_pop, double &intensity) const;
  
    iteration_scheme_line_overlap(const std::string& data_path, int verbosity = 1);
    ~iteration_scheme_line_overlap();
};

// i2 - nb of double overlaps, i3 - nb of triple overlaps
void get_nb_overlap_lines(const energy_diagram* diagram, const einstein_coeff* einst_coeff, double vel_width, int& i2, int& i3);
