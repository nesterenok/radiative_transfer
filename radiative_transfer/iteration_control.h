#pragma once

#include <list>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cstring>

#include "constants.h"
#include "utils.h"
#include "linear_algebra.h"
#include "spectroscopy.h"
#include "coll_rates.h"

void boundary_layer_populations(double *, const energy_diagram*, const einstein_coeff*, const collisional_transitions *, 
	double temp_n, double temp_el, double el_conc, double h_conc, double ph2_conc, double oh2_conc, double he_conc);

// The dynamic array for floating point values;
struct dynamic_array
{
	int dim;
	double *arr;
	
    dynamic_array & operator = (const dynamic_array& obj);
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
	// nb_prev_steps - nb of vectors used;
	// accel_period - nb of iterations between acceleration steps;
	int		nb_prev_steps, accel_start, accel_period, iter_nb, nb_after_accel;
	double	best_eq_error;

	dynamic_array				*opt_level_pop;
	std::list<dynamic_array>	prev_level_pop, residual_list;
	iteration_step				*it_scheme;
	
    // determines new level populations n_{i+1} based on n_i, returns new populations and the parameter is_accelerated
	void next_step(dynamic_array &lev_pop, bool & is_accelerated); 
	void accel_step(double *);
	
public:
	int max_iter_nb, dim;
    // eq_error - maximal abs value in the right hand side of the equation system for level populations, based on n_i (old)
    // pop_error - the maximal abs change in the level populations at the iterative step, max abs(n_{i+1} - n_i) 
    // rel_error - the same, but the relative change,
	double eq_error, pop_error, rel_error;
	
    // this function calculates the populations given the initial populations;
    // the dimension of the population array must be the same as in iteration_step class,
    // min_error is a limit value of the maximal relative population increment between iterations (rel_error),
	bool calculate_populations(double *populations, int max_iter_nb, double min_error, bool acceleration, int verbosity = 1);

    // initialization of parameters: acceleration start, period and nb of vectors used in accelerated step;
	void set_accel_parameters(int start, int period, int nb);

	iteration_control(iteration_step *);
};

template <class iteration_step> 
iteration_control<iteration_step>::iteration_control(iteration_step *its) 
	: it_scheme(its), opt_level_pop(0), nb_prev_steps(5), accel_period(5), accel_start(40), acceleration(false),
	max_iter_nb(0), iter_nb(0), nb_after_accel(0), dim(0), best_eq_error(0.), eq_error(0.), pop_error(0.), rel_error(0.)
{;}

template <class iteration_step> 
void iteration_control<iteration_step>::set_accel_parameters(int start, int period, int nb_p)
{
    acceleration = true;
	accel_start = start;
	accel_period = period;
	nb_prev_steps = nb_p;
}

template <class iteration_step>
void iteration_control<iteration_step>::next_step(dynamic_array & pop_old, bool &is_accelerated)
{
	int i;
	double a;
	dynamic_array residual(dim), pop_new(dim);

	is_accelerated = false;
    prev_level_pop.push_front(pop_old);

	if (acceleration && (iter_nb == accel_start || nb_after_accel == accel_period)) {
		accel_step(pop_old.arr);
		is_accelerated = true;
		nb_after_accel = 0;
	}
	
    // is done in all cases - accelerated step or not, eq_error corresponds to old populations,
	it_scheme->calc_new_pop(pop_old.arr, pop_new.arr, eq_error);	
	if (acceleration && iter_nb >= accel_start)
		nb_after_accel++;
	
    // the level populations are saved that corresponds to a smaller equation error,
	if (eq_error < best_eq_error) {
		best_eq_error = eq_error;
		memcpy(opt_level_pop->arr, pop_old.arr, dim*sizeof(double));
	}
	
	pop_error = rel_error = 0.;
	for (i = 0; i < dim; i++) {
		residual.arr[i] = pop_new.arr[i] - pop_old.arr[i]; // pop_old may be changed in acceleration step,	
		if ((a = fabs(residual.arr[i])) > pop_error) 
            pop_error = a;
		
        if ((a = fabs(residual.arr[i]/(pop_old.arr[i] + 1.e-99))) > rel_error) 
            rel_error = a; 
	}
	residual_list.push_front(residual);
	
	if ((int) residual_list.size() > nb_prev_steps + 1) 
        residual_list.pop_back();
	
    if ((int) prev_level_pop.size() > nb_prev_steps + 1) 
        prev_level_pop.pop_back();

	if (iter_nb < max_iter_nb-1) {
        pop_old = pop_new;
	}
	else {
		// if the population series are not converged, the populations are chosen that have minimum equation error;
        pop_old = *opt_level_pop;
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
	for (i = 0; i < nb_param; i++) {
		i_p++;
		j_p = residual_list.begin();

		for (j = 0; j < nb_param; j++)
		{
			j_p++;
			for (k = 0; k < dim; k++)
			{
				w = prev_level_pop.front().arr[k] + 1.e-99;
                min_arr[i][j] += (residual_list.front().arr[k] - i_p->arr[k])
                    * (residual_list.front().arr[k] - j_p->arr[k])/(w * w);
			}
		}
		for (k = 0; k < dim; k++)
		{
			w = prev_level_pop.front().arr[k] + 1.e-99;
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

// the dimension of the given array has to coincide with the parameter dim of the class;
template <class iteration_step>
bool iteration_control<iteration_step>::calculate_populations(double *pop_arr, int max_nb, double min_error, bool accel, int verbosity)
{
	bool is_found, is_accelerated;
	time_t c_begin = time(NULL);

	acceleration = accel;
	max_iter_nb = max_nb;
	iter_nb = nb_after_accel = 0;
	
	best_eq_error = 1.;
	eq_error = pop_error = rel_error = 0.;

// the dimension of the population arrays is initialized;
	dim = it_scheme->get_vector_dim();
	
    opt_level_pop = new dynamic_array(dim);
	dynamic_array pop_dyn_arr(dim, pop_arr);
	
	std::cout.precision(2);
	std::cout << std::scientific;
	if (verbosity) {
		std::cout << "Level populations calculation:" << std::endl 
            << std::left << std::setw(6) << "Nb" << std::setw(11) << "eq. error " << std::setw(11) << "pop. delta " << std::setw(11) << "rel. pop. delta " << std::endl;
	}
	do {
		next_step(pop_dyn_arr, is_accelerated);
        is_found = (rel_error < min_error);
		
        if (verbosity) {
			if (is_accelerated) 
                std::cout << "acceleration step" << std::endl;
            std::cout << std::left << std::setw(6) << iter_nb << std::setw(11) << eq_error << std::setw(11) << pop_error << std::setw(11) << rel_error << std::endl;
		}
	}
	while (iter_nb < max_iter_nb && !is_found);

    memcpy(pop_arr, pop_dyn_arr.arr, dim*sizeof(double));

	prev_level_pop.clear(); 
	residual_list.clear();
	delete opt_level_pop;
	
	if (verbosity) 
		std::cout << "Time elapsed:	" << (int) (time(NULL) - c_begin) << std::endl;
	return is_found;
}
