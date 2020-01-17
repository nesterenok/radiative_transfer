#pragma once

#include <list>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <cfloat>

#include "constants.h"
#include "utils.h"
#include "linear_algebra.h"
#include "spectroscopy.h"
#include "coll_rates.h"

void boundary_layer_populations(double*, const energy_diagram*, const einstein_coeff*, const collisional_transitions *, 
	double temp_n, double temp_el, double el_conc, double h_conc, double ph2_conc, double oh2_conc, double he_conc);
//void init_slab_popul(const energy_diagram *, const einstein_coeff *, collisional_transitions *, double *populations, int verbosity=1);

// The dynamic array for floating point values;
struct dynamic_array
{
	int dim;
	double *arr;
	
	dynamic_array(int);
	dynamic_array(int, double *);
	~dynamic_array();

	dynamic_array(const dynamic_array &d_a);
};

// The template class that manage the iteration process:
template <class iteration_step> class iteration_control
{
private:
	bool	acceleration;
	// Parameters of the acceleration step: number of vectors used; the minimum iteration number for the accelerated step;
	// the number of iterations between acceleration steps;
	int		nb_prev_steps, default_accel_start, accel_start, accel_period, iter_nb, nb_after_accel;
	double	best_eq_error;

	dynamic_array				*opt_level_pop;
	std::list<dynamic_array>	prev_level_pop, residual_list;
	iteration_step				*it_scheme;
	
    // determines new level populations and the boolean parameter is_accelerated
	void next_step(double *lev_pop, bool & is_accelerated); 
	void accel_step(double *);
	
public:
	int max_iter_nb, dim;
    // eq_error - maximal value in the right hand side of the equation system for level populations
    // pop_error - the maximal change in the level populations at the iterative step, 
    // rel_error - the same, but the relative change,
	double eq_error, pop_error, rel_error;
	
    // this function calculates the populations given the initial populations;
    // the dimension of the population array must be the same as in iteration_step class,
	bool calculate_populations(double *populations, int max_nb, double min_error, bool acceleration = true, int verbosity = 1);

    // initialization of parameters: acceleration start, period and nb of vectors used in accelerated step;
	void set_accel_parameters(int start, int period, int nb);

	iteration_control(iteration_step *);
};

template <class iteration_step> 
iteration_control<iteration_step>::iteration_control(iteration_step *its) 
	: it_scheme(its), opt_level_pop(0), nb_prev_steps(5), accel_period(5), default_accel_start(40), accel_start(0), acceleration(true), 
	max_iter_nb(0), iter_nb(0), nb_after_accel(0), dim(0), best_eq_error(0.), eq_error(0.), pop_error(0.), rel_error(0.)
{;}

template <class iteration_step> 
void iteration_control<iteration_step>::set_accel_parameters(int start, int period, int nb_p)
{
	default_accel_start = accel_start = start;
	accel_period = period;
	nb_prev_steps = nb_p;
}

template <class iteration_step>
void iteration_control<iteration_step>::next_step(double* level_pop, bool &is_accelerated)
{
	int i;
	double a;
	dynamic_array residual(dim), populations_1(dim), populations_2(dim);

	// if the error is small try an acceleration at this step;
	if (rel_error < 3.e-4 && iter_nb > nb_prev_steps+1 && iter_nb < accel_start)
		accel_start = iter_nb;

	is_accelerated = false;
	if (acceleration && (iter_nb == accel_start || nb_after_accel == accel_period))
	{
		accel_step(populations_1.arr);
		is_accelerated = true;
		nb_after_accel = 0;
	}
	else memcpy(populations_1.arr, prev_level_pop.front().arr, dim*sizeof(double));
	
    // is done in all cases - accelerated step or not
	it_scheme->calc_new_pop(populations_1.arr, populations_2.arr, eq_error);	
	if (acceleration && iter_nb >= accel_start)
		nb_after_accel++;
	
// the level populations are saved that correspond to a smaller equation error,
	if (eq_error < best_eq_error) {
		best_eq_error = eq_error;
		memcpy(opt_level_pop->arr, populations_1.arr, dim*sizeof(double));
	}
	
	pop_error = rel_error = 0.;
	for (i = 0; i < dim; i++)
	{
		residual.arr[i] = populations_2.arr[i] - populations_1.arr[i];	
		if ((a = fabs(residual.arr[i])) > pop_error) pop_error = a;
		if ((a = fabs(residual.arr[i]/(populations_1.arr[i] + DBL_EPSILON))) > rel_error) rel_error = a; 
	}

	residual_list.push_front(residual);
	prev_level_pop.push_front(populations_2);
	
	if ((int) residual_list.size() > nb_prev_steps + 1) residual_list.pop_back();
	if ((int) prev_level_pop.size() > nb_prev_steps + 1) prev_level_pop.pop_back();

	if (iter_nb < max_iter_nb-1) {
		memcpy(level_pop, populations_2.arr, dim*sizeof(double));
	}
	else {
		// If the population series are not converged, the populations are chosen that have minimum equation error;
		memcpy(level_pop, opt_level_pop->arr, dim*sizeof(double));
		eq_error = best_eq_error;
	}
	iter_nb++;
}

template <class iteration_step>
void iteration_control<iteration_step>::accel_step(double* accel_pop)
{
	int i, j, k, nb_param;
	double w, sum, *min_b, **min_arr;	
	
	nb_param = nb_prev_steps-1; // the number of unknowns is less by 1 than the number of population vectors used;
	min_b = new double [nb_param];
	min_arr = alloc_2d_array<double>(nb_param, nb_param);

	std::list<dynamic_array>::iterator i_p, j_p;

	memset(*min_arr, 0, nb_param*nb_param*sizeof(double));
	memset(min_b, 0, nb_param*sizeof(double));

	i_p = residual_list.begin();
	for (i = 0; i < nb_param; i++)
	{
		i_p++;
		j_p = residual_list.begin();

		for (j = 0; j < nb_param; j++)
		{
			j_p++;
			for (k = 0; k < dim; k++)
			{
				w = prev_level_pop.front().arr[k] + DBL_EPSILON;
				min_arr[i][j] += (residual_list.front().arr[k] - i_p->arr[k]) 
					*(residual_list.front().arr[k] - j_p->arr[k])/(w*w);
			}
		}
		for (k = 0; k < dim; k++)
		{
			w = prev_level_pop.front().arr[k] + DBL_EPSILON;
            min_b[i] += (residual_list.front().arr[k] - i_p->arr[k]) * residual_list.front().arr[k]/ (w * w);
		}
	}
	
	lu_matrix_solve(min_arr, min_b, nb_param);
	
	sum = 0.;
	for (i = 0; i < nb_param; i++) {
		sum += min_b[i];
	}
	for (k = 0; k < dim; k++) {	
		i_p = prev_level_pop.begin();
		accel_pop[k] = (1. - sum)* i_p->arr[k];
		
		for (i = 0; i < nb_param; i++) {	
			i_p++;
			accel_pop[k] += min_b[i]* i_p->arr[k];
		}
	}
	delete [] min_b;
	free_2d_array<double>(min_arr);
}

// min_error is a limit value of the maximal relative population increment between iterations (rel_error),
// the dimension of the given array has to coincide with the parameter dim of the class;
// During the calculations, all physical parameters have not to be changed;
template <class iteration_step>
bool iteration_control<iteration_step>::calculate_populations(double *populations, int max_nb, double min_error, bool accel, int verbosity)
{
	bool	is_found, is_accelerated;
	time_t	c_begin = time(NULL);

	acceleration = accel;
	max_iter_nb = max_nb;
	iter_nb = nb_after_accel = 0;
	accel_start = default_accel_start;
	
	best_eq_error = 1.;
	eq_error = pop_error = rel_error = 0.;

// the dimension of the population arrays is initialized;
	dim = it_scheme->get_vector_dim();
	opt_level_pop = new dynamic_array(dim);

// Initial populations are added to the list with level population arrays; 
	dynamic_array init_pop(dim, populations);
	prev_level_pop.push_front(init_pop);

	std::cout.precision(3);
	std::cout << std::scientific;
	if (verbosity) {
		std::cout << "Level populations calculation:" << std::endl << std::left
			<< std::setw(5) << "Nb" << std::setw(15) << "eq. error" << std::setw(15) << "pop. delta" << std::setw(15) 
			<< "rel. pop. delta" << std::endl;
	}
	do {
		next_step(populations, is_accelerated);
		if (verbosity) 
		{
			if (is_accelerated) std::cout << "acceleration step" << std::endl;
			std::cout << std::left << std::setw(5) << iter_nb << std::setw(15) << eq_error << std::setw(15) << pop_error 
				<< std::setw(15) << rel_error << std::endl;
		}
		is_found = (rel_error < min_error);
	}
	while (iter_nb < max_iter_nb && !is_found);

// After calculations, all auxiliary data are deleted;	
	prev_level_pop.clear(); 
	residual_list.clear();
	delete opt_level_pop;
	
	if (verbosity) 
		std::cout << "Time elapsed:	" << (int) (time(NULL) - c_begin) << std::endl;
	return is_found;
}
