
#include "transition_list.h"
using namespace std;

// The list of CH3OH transitions for which the parameters are calculated, class I masers,
// the compilation of observed maser transitions is given by Ladeyschikov et al., Astronomical J. 158, 233 (2019), 
void ch3oh_classI_trans_list(const energy_diagram* diagram, list<transition>& trans_list)
{
    int	up, low;
    transition* trans;

    // A methanol species
    if (rounding(2. * diagram->mol.spin) == 3)
    {
        // K < 0 means A- level, 
        up = diagram->get_nb(0, 10, -1);	//  23.444 GHz; discovery: Voronkov et al. MNRAS 413, 2339, 2011
        low = diagram->get_nb(0, 9, -2);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 7, 0);	//  44.069 GHz, detected
        low = diagram->get_nb(0, 6, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 8, 0);	//  95.169 GHz, detected
        low = diagram->get_nb(0, 7, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 9, 0);	//  146.618 GHz
        low = diagram->get_nb(0, 8, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 10, 0);	//  198.403 GHz
        low = diagram->get_nb(0, 9, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 11, 0);	//  250.507 GHz
        low = diagram->get_nb(0, 10, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }
    }

    // E methanol species:
    else if (rounding(2. * diagram->mol.spin) == 1)
    {
        up = diagram->get_nb(0, 9, -1);	//  9.9362 GHz, detected
        low = diagram->get_nb(0, 8, -2);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        // there is a set of maser transitions J_2->J_1, J = 2,3,4,5,6,7,8,9 (Voronkov et al., MNRAS 373, 411-424, 2006),
        // frequencies from Mekhtiev et al., J. Mol. Spectr. 194, 171-178 (1999),
        up = diagram->get_nb(0, 2, 2);	//  24.934398 GHz
        low = diagram->get_nb(0, 2, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 3, 2);	//  24.928715 GHz
        low = diagram->get_nb(0, 3, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 4, 2);	//  24.93348 GHz
        low = diagram->get_nb(0, 4, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 5, 2);	//  24.959084 GHz
        low = diagram->get_nb(0, 5, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 6, 2);	//  25.018122 GHz
        low = diagram->get_nb(0, 6, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 7, 2);	//  25.124864 GHz
        low = diagram->get_nb(0, 7, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 8, 2);	//  25.294401 GHz
        low = diagram->get_nb(0, 8, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 9, 2);	//  25.541375 GHz
        low = diagram->get_nb(0, 9, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 10, 2);	//  25.878239 GHz
        low = diagram->get_nb(0, 10, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }
        // end of J_2->J_1 series

        up = diagram->get_nb(0, 4, 0);	//  28.316 GHz
        low = diagram->get_nb(0, 3, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 4, -1);	//  36.169 GHz, detected
        low = diagram->get_nb(0, 3, 0);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 5, 0);	//  76.509 GHz
        low = diagram->get_nb(0, 4, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 5, -1);	//  84.521 GHz, detected
        low = diagram->get_nb(0, 4, 0);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 11, -1); //  104.300 GHz, detected
        low = diagram->get_nb(0, 10, -2);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 6, 0); //  124.569 GHz
        low = diagram->get_nb(0, 5, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 6, -1); //  132.890 GHz, detected
        low = diagram->get_nb(0, 5, 0);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 7, -1); //  181.296 GHz, detected
        low = diagram->get_nb(0, 6, 0);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 4, 2); //  218.440 GHz, detected
        low = diagram->get_nb(0, 3, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 8, -1); //  229.758 GHz, detected
        low = diagram->get_nb(0, 7, 0);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 5, 2); //  266.838 GHz, detected
        low = diagram->get_nb(0, 4, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 9, -1); //  278.304 GHz, detected
        low = diagram->get_nb(0, 8, 0);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }
    }
}

// The list of CH3OH transitions for which the parameters are calculated, class II masers
void ch3oh_classII_trans_list(const energy_diagram* diagram, list<transition>& trans_list)
{
    int	up, low;
    transition* trans;

    // A methanol species
    if (rounding(2. * diagram->mol.spin) == 3)
    {
        up = diagram->get_nb(0, 5, 1);	//  6.669 GHz
        low = diagram->get_nb(0, 6, 0);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 9, 2);	//  23.121 GHz
        low = diagram->get_nb(0, 10, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 6, -2);	//  38.293 GHz
        low = diagram->get_nb(0, 5, -3);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 6, 2);	//  38.453 GHz
        low = diagram->get_nb(0, 5, 3);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 4, 1);	//  57.033 GHz
        low = diagram->get_nb(0, 5, 0);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 7, -2);	//  86.616 GHz
        low = diagram->get_nb(0, 6, -3);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 7, 2);	//  86.903 GHz
        low = diagram->get_nb(0, 6, 3);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 3, 1);	//  107.013 GHz
        low = diagram->get_nb(0, 4, 0);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }
    }
    // E methanol species:
    else if (rounding(2. * diagram->mol.spin) == 1)
    {
        up = diagram->get_nb(0, 2, 0);	// 12.179 GHz
        low = diagram->get_nb(0, 3, -1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 2, 1);	// 19.967 GHz; Krishnan et al., MNRAS 433, 3346, 2013
        low = diagram->get_nb(0, 3, 0);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 7, -2);	// 37.704 GHz
        low = diagram->get_nb(0, 8, -1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 6, -2);	// 85.568 GHz
        low = diagram->get_nb(0, 7, -1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 0, 0);	//  108.894 GHz
        low = diagram->get_nb(0, 1, -1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        // There is a set of J,0->J,-1 transitions:
        up = diagram->get_nb(0, 5, 0);	//  157.179 GHz
        low = diagram->get_nb(0, 5, -1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        // Torsionally excited:
        up = diagram->get_nb(1, 2, 0);	//  44.956 GHz; Voronkov et al., A&A 387, 310, 2002
        low = diagram->get_nb(1, 3, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }
    }
}


void nh3_trans_list(const energy_diagram* diagram, std::list<transition>& trans_list)
{
    int up, low;
    transition* trans;

    if (rounding(2. * diagram->mol.spin) == 3) // ortho 
    {
        up = diagram->get_nb(1, 0, 3, 3);	//  23.870 GHz = 0.79622 cm-1
        low = diagram->get_nb(0, 0, 3, 3);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }
    }

    if (rounding(2. * diagram->mol.spin) == 1) // para
    {
        up = diagram->get_nb(1, 0, 2, 2);	//  23.723 GHz = 0.7913 cm-1
        low = diagram->get_nb(0, 0, 2, 2);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }
    }
}

void oh_trans_list(const energy_diagram* diagram, std::list<transition>& trans_list)
{
    int up, low;
    transition* trans;

    // parity, v, j, omega, F
    up = diagram->get_nb(1, 0, 1.5, 1.5, 2.);	//  1721 MHz
    low = diagram->get_nb(-1, 0, 1.5, 1.5, 1.);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(1, 0, 1.5, 1.5, 1.);	//  1612 MHz
    low = diagram->get_nb(-1, 0, 1.5, 1.5, 2.);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(-1, 0, 2.5, 1.5, 3.);	//  6050 MHz
    low = diagram->get_nb(1, 0, 2.5, 1.5, 2.);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(-1, 0, 0.5, 0.5, 1.);	//  4765 MHz
    low = diagram->get_nb(1, 0, 0.5, 0.5, 0.);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }
}

void h2co_trans_list(const energy_diagram *diagram, std::list<transition>& trans_list)
{
    int up, low;
    transition* trans;

    // v, j, ka - kc
    up = diagram->get_nb(0, 1, 1);	//  4.83 GHz
    low = diagram->get_nb(0, 1, 0);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(0, 2, 0);	//  14.5 GHz
    low = diagram->get_nb(0, 2, -1);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }
}

void para_h2o_trans_list(const energy_diagram *diagram, std::list<transition>& trans_list)
{
    int up, low;
    transition* trans;

    up = diagram->get_nb(0, 3, -2);	// 183.3 GHz
    low = diagram->get_nb(0, 2, 2);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(0, 5, -4); // 325.2 GHz
    low = diagram->get_nb(0, 4, 0);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    // Vibrationally excited
    up = diagram->get_nb(1, 4, 4);	// 96.26 GHz
    low = diagram->get_nb(1, 5, 0);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }
}

void ortho_h2o_trans_list(const energy_diagram* diagram, std::list<transition>& trans_list)
{
    int up, low;
    transition* trans;

    up = diagram->get_nb(0, 6, -5);	// 22.235 GHz
    low = diagram->get_nb(0, 5, -1);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(0, 10, -7); // 321.2 GHz
    low = diagram->get_nb(0, 9, -3);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(0, 17, -9); // 354.8 GHz
    low = diagram->get_nb(0, 16, -3);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(0, 4, -3);	// 380.2 GHz
    low = diagram->get_nb(0, 3, 1);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(0, 6, 1);	// 439.15 GHz
    low = diagram->get_nb(0, 5, 5);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(0, 7, 3);	// 443.0 GHz
    low = diagram->get_nb(0, 6, 5);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(0, 4, -1);	// 448.0 GHz
    low = diagram->get_nb(0, 3, 3);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(0, 5, 1);	// 620.7 GHz
    low = diagram->get_nb(0, 4, 3);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }
}

void ortho_h2o_vibr_trans_list(const energy_diagram *diagram, std::list<transition> & trans_list)
{
    int up, low;
    transition* trans;

    // Vibrationally excited
    up = diagram->get_nb(1, 4, -1);	// 12.01 GHz
    low = diagram->get_nb(1, 3, 3);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(1, 4, -3);	// 67.80 GHz
    low = diagram->get_nb(1, 3, 1);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(1, 5, 5);	// 232.7 GHz
    low = diagram->get_nb(1, 6, 1);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(1, 6, 5);	// 293.7 GHz
    low = diagram->get_nb(1, 7, 3);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(1, 5, -1);	// 336.2 GHz
    low = diagram->get_nb(1, 6, -5);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(3, 1, 1);	// 540.7 GHz
    low = diagram->get_nb(3, 1, -1);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(1, 1, 1);	// 658.0 GHz
    low = diagram->get_nb(1, 1, -1);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }

    up = diagram->get_nb(2, 1, 1);	// 793.6 GHz
    low = diagram->get_nb(2, 1, -1);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }
}
