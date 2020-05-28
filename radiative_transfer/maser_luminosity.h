#pragma once
#include "iteration_lvg.h"
#include "transition_data.h"

// Calculates the limiting luminosity, pump efficiency, loss rate,
// limiting luminosity - the rate of maser photon production in unit volume [cm-3 s-1],
// pump efficiency = inversion/(n_i/g_i + n_j/g_j)
// if there is no population inversion, the parameters equals to 1.e-99;
// the answer is saved to the object of transition_data_container,
void lim_luminosity_lvg(iteration_scheme_lvg* calc_scheme, transition_data_container* container, const cloud_data* cloud, 
    const energy_diagram* diagram, const einstein_coeff* einst_coeff, collisional_transitions* c_trans, double* level_pop);
