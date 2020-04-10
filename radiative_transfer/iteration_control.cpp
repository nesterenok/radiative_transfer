
#include <cmath>
#include <float.h>
#include <memory.h>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "iteration_control.h"
#include "utils.h"
#include "linear_algebra.h"

using namespace std;

dynamic_array::dynamic_array(int d)
{
	dim = d;
	arr = new double [dim];
	memset(arr, 0, dim*sizeof(double));
}

dynamic_array::dynamic_array(int d, double *arr2)
{
	dim = d;
	arr = new double [dim];
	memcpy(arr, arr2, dim*sizeof(double));
}

dynamic_array::~dynamic_array() {
	delete [] arr;
}

dynamic_array::dynamic_array(const dynamic_array &obj)
{
	dim = obj.dim;
	arr = new double [dim];
	memcpy(arr, obj.arr, dim*sizeof(double));
}

dynamic_array& dynamic_array::operator = (const dynamic_array& obj)
{
    dim = obj.dim;
    delete[] arr;
    arr = new double[dim];
    
    memcpy(arr, obj.arr, dim * sizeof(double));  
    return *this;
}

// The physical parameters in collisional_transitions object have to be initialized before the function call; 
void boundary_layer_populations(double *arr, const energy_diagram *diagram, const einstein_coeff *e_coeff, const collisional_transitions *coll_trans,
	double temp_n, double temp_el, double el_conc, double h_conc, double ph2_conc, double oh2_conc, double he_conc)
{
	int	i, j, nb_mol_lev;
	int *indices(0);
	
	double down, up;
	double *coll_partn_conc(0), **matrix;

	coll_trans->set_gas_param(temp_n, temp_el, he_conc, ph2_conc, oh2_conc, h_conc, el_conc, coll_partn_conc, indices);

	nb_mol_lev = diagram->nb_lev;
	matrix = alloc_2d_array<double>(nb_mol_lev, nb_mol_lev);

	memset(*matrix, 0, nb_mol_lev*nb_mol_lev*sizeof(double));
	memset(arr, 0, nb_mol_lev*sizeof(double));
	
	for (i = 1; i < nb_mol_lev; i++) {
		for (j = 0; j < i; j++)
		{
			coll_trans->get_rate_neutrals(diagram->lev_array[i], diagram->lev_array[j], down, up, temp_n, coll_partn_conc, indices);

			matrix[j][i] = 0.5*e_coeff->arr[i][j] + down;	// i->j
			matrix[i][i] -= 0.5*e_coeff->arr[i][j] + down;
		
			matrix[i][j] = up;	// j->i
			matrix[j][j] -= up;
		}
	}

	for (i = 0; i < nb_mol_lev; i++) {
		matrix[0][i] = 1.;
	}
	arr[0] = 1.;
	
	lu_matrix_solve(matrix, arr, nb_mol_lev);
	free_2d_array<double>(matrix);
    delete[] indices;
    delete[] coll_partn_conc;
}

/*
//void init_slab_popul(const energy_diagram *, const einstein_coeff *, collisional_transitions *, double *populations, int verbosity=1);

// The function calculates the initial populations. The cloud is approximated by a single layer. 
// The iteration technique is implemented. The initial guess for populations is Boltzmann or "boundary layer" distribution.
void init_slab_popul(const slab_cloud *sl_cloud, const energy_diagram *diagram, const einstein_coeff *einst_coeff, 
	collisional_transitions *coll_trans, double* populations, int verbosity)
{
	bool	is_found, acceleration;
	int		lay, nb_cloud_lay, nb_mol_lev, iter_max_nb;
	double	min_error, cloud_height;
	double	*bound_lay_pop, *aux_pop_arr;
	
	const cloud_layer *layer;
	
	iter_max_nb = 30;
	min_error = 1.e-2;
	acceleration = true;
	
	nb_mol_lev = diagram->nb_lev;
	nb_cloud_lay = sl_cloud->nb_lay;
	cloud_height = sl_cloud->height;

	bound_lay_pop = new double [nb_mol_lev];
	aux_pop_arr = new double [nb_mol_lev];
	if (verbosity) cout << "Initialization of populations:" << endl;
	
	// Initialization of the "one-layer" iterative scheme, the external radiation is not included!!
	escape_probability_fixed_depth *it_scheme = new escape_probability_fixed_depth();
	
	it_scheme->set_dust_model(sl_cloud->dust);
	it_scheme->init_molecule_data(diagram, einst_coeff, coll_trans);

	iteration_control<escape_probability_fixed_depth> it_control(it_scheme);
	
	// init parameters for acceleration step in iteration series: start, period and nb of vectors used in acceleration step;
	it_control.set_accel_parameters(6, 4, 4);
	
	for (lay = 0; lay < nb_cloud_lay; lay++)
	{
		layer = &sl_cloud->lay_array[lay];
		// The parameters of the collisional_transitions object are initialized here:
		coll_trans->set_phys_param(layer->gas_temp, layer->he_conc, layer->ph2_conc, layer->oh2_conc, layer->h_conc, layer->el_conc);

		if (lay == 0) {
			boundary_layer_populations(bound_lay_pop, diagram, einst_coeff, coll_trans);
			memcpy(aux_pop_arr, bound_lay_pop, nb_mol_lev*sizeof(double));
		}
		// The parameters of the iterative scheme are initialized:
		it_scheme->set_depth(layer->z_m);
		it_scheme->set_parameters(cloud_height, layer->gas_temp, layer->dust_temp, layer->el_conc, layer->h_conc, 
			layer->ph2_conc, layer->oh2_conc, layer->he_conc, layer->mol_conc, layer->vel_turb);
		
		is_found = it_control.calculate_populations(aux_pop_arr, iter_max_nb, min_error, acceleration, false);	
		if (is_found) {
			memcpy(populations + lay*nb_mol_lev, aux_pop_arr, nb_mol_lev*sizeof(double));
		}
		else {
			boltzmann_populations(layer->gas_temp, aux_pop_arr, diagram);
			is_found = it_control.calculate_populations(aux_pop_arr, iter_max_nb, min_error, acceleration, false);
			
			if (is_found) {
				memcpy(populations + lay*nb_mol_lev, aux_pop_arr, nb_mol_lev*sizeof(double));
			}
			else {
				if (verbosity) cout << left << setw(5) << lay << "n";
				if (lay > 0) memcpy(populations + lay*nb_mol_lev, populations + (lay-1)*nb_mol_lev, nb_mol_lev*sizeof(double));		
				else memcpy(populations, bound_lay_pop, nb_mol_lev*sizeof(double));
			}
		}
	
		if (is_found && verbosity) cout << left << setw(5) << lay;
		if (verbosity && (lay+1)%15 == 0) cout << endl; 
	}
	if (verbosity) cout << endl;

	delete [] aux_pop_arr;
	delete [] bound_lay_pop;
	
	delete it_scheme;
}*/
