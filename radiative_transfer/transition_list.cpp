
#include "transition_list.h"
using namespace std;

// The list of CH3OH transitions for which the parameters are calculated, class I masers:
void ch3oh_classI_trans_list(const energy_diagram* diagram, list<transition>& trans_list)
{
    int	up, low;
    transition* trans;

    // A methanol species
    if (rounding(2. * diagram->mol.spin) == 3)
    {
        up = diagram->get_nb(0, 10, -1);	//  23.445 GHz; discovery: Voronkov et al. MNRAS 413, 2339, 2011
        low = diagram->get_nb(0, 9, -2);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 7, 0);	//  44.069 GHz
        low = diagram->get_nb(0, 6, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 8, 0);	//  95.169 GHz
        low = diagram->get_nb(0, 7, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }
    }

    // E methanol species:
    else if (rounding(2. * diagram->mol.spin) == 1)
    {
        up = diagram->get_nb(0, 9, -1);	//  9.9362 GHz
        low = diagram->get_nb(0, 8, -2);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        // There is a series of maser transitions J_2->J_1, J = 2,3,4,5,6,7,8,9
        up = diagram->get_nb(0, 5, 2);	//  24.959 GHz
        low = diagram->get_nb(0, 5, 1);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 4, -1);	//  36.169 GHz
        low = diagram->get_nb(0, 3, 0);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 5, -1);	//  84.521 GHz
        low = diagram->get_nb(0, 4, 0);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }

        up = diagram->get_nb(0, 11, -1);	//  104.30 GHz
        low = diagram->get_nb(0, 10, -2);

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

    if (rounding(2. * diagram->mol.spin) == 3)
    {
        up = diagram->get_nb(1, 0, 3, 3);	//  23.870 GHz = 0.79622 cm-1
        low = diagram->get_nb(0, 0, 3, 3);

        if (up != -1 && low != -1) {
            trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
            trans_list.push_back(*trans);
            delete trans;
        }
    }

    if (rounding(2. * diagram->mol.spin) == 1)
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

/*    up = diagram->get_nb(1, 0, 1.5, 1.5, 1.);	//  1612 MHz
    low = diagram->get_nb(-1, 0, 1.5, 1.5, 2.);

    if (up != -1 && low != -1) {
        trans = new transition(diagram->lev_array[low], diagram->lev_array[up]);
        trans_list.push_back(*trans);
        delete trans;
    }*/

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
