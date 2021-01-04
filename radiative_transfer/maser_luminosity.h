#pragma once
#include "iteration_lvg.h"
#include "transition_data.h"

// Calculates the limiting luminosity, pump efficiency, loss rate, emission measure,
// limiting luminosity = the rate of maser photon production in unit volume [cm-3 s-1], (the definition is different than in Leurini et al. 2016)
// pump efficiency = inversion/(n_i/g_i + n_j/g_j)
// pump rate = 0.5*(pump_rate_up/g_up + pump_rate_low/g_low)/(n_H2 *n_mol), pump rate is per level, pump_rate = population density *loss_rate  [cm-3 s-1]
// loss_rate = (up_loss_rate* g_up + low_loss_rate* g_low)/(g_up + g_low)  [s-1]
// emission measure = n_H2 *n_mol/(dv/dz)
// if there is no population inversion, the parameters equals to 1.e-99; the answer is saved to the object of transition_data_container,
void lim_luminosity_lvg(iteration_scheme_lvg* calc_scheme, transition_data_container* container, const cloud_data* cloud, 
    const energy_diagram* diagram, const einstein_coeff* einst_coeff, collisional_transitions* c_trans, double* level_pop);
