
#include <cstring>
#include "maser_luminosity.h"
using namespace std;

// It is assumed that: the molecular data, radiation field, dust model are initialized in calc_scheme,
void lim_luminosity_lvg(iteration_scheme_lvg* calc_scheme, transition_data_container* container, const cloud_data* cloud, const energy_diagram* diagram,
    const einstein_coeff* einst_coeff, collisional_transitions* coll_trans, double* level_pop)
{
    int nb_cloud_lay, nb_mol_lev, lay, i, j, k, shift;
    double inversion, intensity, loss_rate;

    int* indices(0);
    double* up_loss_rate, * low_loss_rate, * coll_partn_conc(0);

    energy_level l_low, l_up;
    cloud_layer clayer;
    std::list<transition_data>::iterator it;

    nb_cloud_lay = cloud->nb_lay;
    nb_mol_lev = diagram->nb_lev;

    up_loss_rate = new double[nb_cloud_lay];
    memset(up_loss_rate, 0, nb_cloud_lay * sizeof(double));

    low_loss_rate = new double[nb_cloud_lay];
    memset(low_loss_rate, 0, nb_cloud_lay * sizeof(double));

    it = container->data.begin();
    while (it != container->data.end()) {
        l_low = it->trans->low_lev;
        l_up = it->trans->up_lev;
        it->lum = 0.;

        for (lay = 0; lay < nb_cloud_lay; lay++) {
            shift = nb_mol_lev * lay;
            clayer = cloud->lay_array[lay];

            calc_scheme->set_vel_grad(clayer.velg_n);
            calc_scheme->set_dust_parameters(clayer.dust_grain_conc, clayer.dust_grain_temp);
            calc_scheme->set_parameters(clayer.temp_n, clayer.temp_el, clayer.el_conc, clayer.h_conc, clayer.ph2_conc, clayer.oh2_conc,
                clayer.he_conc, clayer.mol_conc, clayer.vel_turb);

            for (k = 0; k <= 1; k++) {
                if (k == 0)
                    j = l_low.nb;
                else j = l_up.nb;

                loss_rate = 0.;
                for (i = 0; i < nb_mol_lev; i++) {
                    if (einst_coeff->arr[i][j] != 0 && i != l_low.nb && i != l_up.nb)
                    { // there is no population inversion between i and j levels, a_ji = a_ij *g_i/g_j
                        if ((i < j) && (level_pop[shift + i] * einst_coeff->arr[i][j] > level_pop[shift + j] * einst_coeff->arr[j][i])) {
                            calc_scheme->intensity_calc(j, i, level_pop, intensity);
                            loss_rate += einst_coeff->arr[j][i] * (1. + intensity);
                        }
                        else if ((i > j) && (level_pop[shift + i] * einst_coeff->arr[i][j] < level_pop[shift + j] * einst_coeff->arr[j][i])) { 
                            calc_scheme->intensity_calc(i, j, level_pop, intensity);
                            loss_rate += einst_coeff->arr[j][i] * intensity;
                        }
                    }
                }
                if (k == 0)
                    low_loss_rate[lay] = loss_rate;
                else up_loss_rate[lay] = loss_rate;
            }

            coll_trans->set_gas_param(clayer.temp_n, clayer.temp_el, clayer.he_conc, clayer.ph2_conc, clayer.oh2_conc, clayer.h_conc, clayer.el_conc,
                coll_partn_conc, indices);

            for (i = 0; i < nb_mol_lev; i++) {
                if (i != l_low.nb)
                    low_loss_rate[lay] += coll_trans->get_rate_neutrals(l_low, diagram->lev_array[i], clayer.temp_n, coll_partn_conc, indices);

                if (i != l_up.nb)
                    up_loss_rate[lay] += coll_trans->get_rate_neutrals(l_up, diagram->lev_array[i], clayer.temp_n, coll_partn_conc, indices);
            }

            // emission measure,
            it->emiss_coeff_arr[lay] = (clayer.ph2_conc + clayer.oh2_conc) * clayer.mol_conc / clayer.velg_n;

            // if inversion < 0 than luminosity = 1.e-99;
            inversion = level_pop[shift + l_up.nb] / l_up.g - level_pop[shift + l_low.nb] / l_low.g;
            if (inversion > 0.) {
                it->lum_arr[lay] = inversion / (1. / (up_loss_rate[lay] * l_up.g) + 1. / (low_loss_rate[lay] * l_low.g)) * clayer.mol_conc;
               
                it->pump_eff_arr[lay] = inversion / (level_pop[shift + l_up.nb] / l_up.g + level_pop[shift + l_low.nb] / l_low.g);
            }
            else {
                it->lum_arr[lay] = it->pump_eff_arr[lay] = 1.e-99;
            }
             
            it->lum += it->lum_arr[lay] * cloud->lay_array[lay].dz;
            it->loss_rate_arr[lay] = (up_loss_rate[lay]* l_up.g + low_loss_rate[lay]* l_low.g)/((double) l_up.g + l_low.g);
            
            it->pump_rate_arr[lay] = 0.5 * (level_pop[shift + l_up.nb] * up_loss_rate[lay] + level_pop[shift + l_low.nb] * low_loss_rate[lay])
                    / (clayer.ph2_conc + clayer.oh2_conc);  // level populations are normalized on molecular concentration
        }      
        it->lum /= cloud->get_height();
        it++;
    }
    delete[] indices;
    delete[] coll_partn_conc;
    delete[] up_loss_rate;
    delete[] low_loss_rate;
}
