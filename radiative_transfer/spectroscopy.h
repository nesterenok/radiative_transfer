#pragma once
#include <string>
#include <vector>
#include "utils.h"

// Max nb of excited vibrational states of H2O is 4:
// the collisional data for rovibrational transitions is given by Faure & Josselin, A&A 492, 257–264 (2008);
#define NB_VIBR_EXCIT_H2O 4

// Nb of torsionally excited states of methanol, for which the radiative coefficients were calculated by Mekhtiev et al., J. Mol. Spectr. 194, 171-178 (1999):
#define NB_VIBR_EXCIT_CH3OH_MEKHTIEV 2
// the maximal value of the angular momentum J for the data by Mekhtiev et al. (1999):
#define NB_ANG_MOM_CH3OH 22

// Nb of vibrationally excited states of CO; 
// collisional data for rovibrational transitions for CO-H collisions is given by Song et al., ApJ 813, 96, (2015);
#define NB_VIBR_EXCIT_CO 5

class energy_diagram;
class einstein_coeff;

// The function calculates thermodynamic equilibrium populations of the molecule;
void boltzmann_populations(double, double*, const energy_diagram*);

// The function calculates the depth in cm, at which optical depth in all lines of specific molecule is less than or equal unity,
// molecule concentration is taken to be 1 cm^{-3};
double calc_depth(const energy_diagram *, const einstein_coeff *, double gas_temp, double vel_disp);

void save_populations(const std::string &output_path, const energy_diagram *diagram, double *arr, int nb_cloud_lay, 
	int nb_mol_lev, bool normalized = true, std::string id = "");

void read_populations(const std::string &file_name, double *arr, int nb_cloud_lay, int nb_mol_lev);

// The class contains the molecule data,
class molecule
{
public:
	int isotop;
	double spin, mass; // the mass must be in [g]
	std::string	name; // the symbol defining the spin isomer must be included in the name
	
	molecule(const std::string& n, int is, double m, double sp = -1.);
};

// The class contains energy level data
class energy_level
{
public:
	// level nb, sym inv (+1,-1), statistical weight (including j and spin), vibrational quantum number, 
	int nb, syminv, g, v;
	// angular momentum, projections of the angular momentum on molecule axes, molecule spin, energy (cm-1), hyperfine splitting quantum number
	double j, k1, k2, spin, energy, hf;
	std::string name;

	// The relation operators are needed to sort by energy; the levels are equal if all quantum numbers coincide,
    // for the ions the variables "name" and "g" define the level 
	bool operator == (const energy_level &obj) const {
		return (v == obj.v && syminv == obj.syminv && rounding(2.*j) == rounding(2.*obj.j) && rounding(2.*k1) == rounding(2.*obj.k1) 
			&& rounding(2.*k2) == rounding(2.*obj.k2) && rounding(2.*hf) == rounding(2.*obj.hf) 
            && name == obj.name && g == obj.g);
	}
	bool operator != (const energy_level &obj) const {return !(*this == obj);}
	bool operator < (const energy_level &obj) const {return (energy < obj.energy && !(*this == obj));}	
	bool operator > (const energy_level &obj) const {return (energy > obj.energy && !(*this == obj));}
	energy_level();
};

// The class contains the energy level table;
class energy_diagram
{
public:
    bool hyperfine_splitting;
	int nb_lev, verbosity;
	const molecule mol;
	std::vector<energy_level> lev_array;

	virtual int get_nb(int v, double j) const {return -1;}
	virtual int get_nb(int v, double j, double k) const {return -1;}
	virtual int get_nb(const std::string, int stat_w) const {return -1;}
	virtual int get_nb(int syminv, int v, double j, double k) const {return -1;}
    virtual int get_nb(int syminv, int v, double j, double k, double hf) const { return -1; }
    
    void report(const std::string & fname);

	energy_diagram(molecule m, int vebosity = 1);
	virtual ~energy_diagram();
};

// H2 molecule (Dabrowski I., Can. J. Phys. 62, p. 1639, 1984); maximal value of rovibrational levels is 318;
class h2_diagram : public energy_diagram
{
public:
	h2_diagram(const std::string &path, molecule m, int &n_l, int verbosity = 1);
	int get_nb(int v, double j) const;
};

// H2O molecule (HITRAN2012 database);
// if nb of vibrationally excited states equals 4, maximal nb of levels is 411, equals 1 - 271 levels, equals 0 - 164 levels;
// either ortho- or para-H2O molecules are considered - molecule spin must be defined,
class h2o_diagram : public energy_diagram
{
public:
	int get_nb(int v, double j, double tau) const;
	int get_vibr_nb(int v1, int v2, int v3) const;
	h2o_diagram(const std::string &path, molecule m, int &n_l, int nb_vibr =NB_VIBR_EXCIT_H2O, int verbosity =1);
};

// CH3OH molecule (Mekhtiev et al., Journal of Molecular Spectroscopy 194, 171-178, 1999);
// the maximum number of levels with vt <= 2 and angular momentum <= 22 is 1455;
// either A or E CH3OH spin isomers is considered - molecule spin must be defined, A type - spin 3/2, E type - spin 1/2;
// A- and A+ differ in the sign of k1,
class ch3oh_diagram : public energy_diagram
{
public:
	int get_nb(int v, double j, double k) const;
	ch3oh_diagram(const std::string &path, molecule m, int &n_l, int nb_vibr =NB_VIBR_EXCIT_CH3OH_MEKHTIEV, int ang_mom_max =NB_ANG_MOM_CH3OH, int verbosity =1);
};

// Energy level table of specific ion. Please, check the names of ion and data file;
// NIST https://www.nist.gov/pml/atomic-spectra-database
class ion_diagram : public energy_diagram
{
public:
	// the name of the level is not unique:
	int get_nb(const std::string, int stat_w) const;
	ion_diagram(const std::string &path, molecule m, int &n_l, int verbosity =1);
};

// CO molecule energy levels. HIRTRAN2008 data.
class co_diagram : public energy_diagram
{
public:
	int get_nb(int v, double j) const;
	co_diagram(const std::string &path, molecule m, int &n_l, int nb_vibr =NB_VIBR_EXCIT_CO, int verbosity =1);
};

// v, j, k1 is omega, syminv is parity, 
class oh_diagram : public energy_diagram
{
public:
	int get_nb(int parity, int v, double j, double omega) const;
	oh_diagram(const std::string &path, molecule m, int &n_l, int verbosity =1);
};

// the hyperfine splitting of OH is taken into account in this class:
// additional quantum number hf - total angular momentum
class oh_hf_diagram : public energy_diagram
{
public:
    // hf is the total angular momentum F = J + S
    int get_nb(int parity, int v, double j, double omega, double hf) const;
    oh_hf_diagram(const std::string& path, molecule m, int& n_l, int verbosity = 1);
};

// NH3, either ortho- or para- molecules are considered - molecule spin must be defined,
class nh3_diagram : public energy_diagram
{
public:
	int get_nb(int syminv, int v, double j, double k) const;
	nh3_diagram(const std::string &path, molecule m, int &n_l, int verbosity =1);
};

// H2CO, either ortho- or para- molecules are considered - molecule spin must be defined,
// some para- levels have equal energy - must be taken into account in the code or fixed in file data, the level with higher k2 must have lower nb, 
class h2co_diagram : public energy_diagram
{
public:
	int get_nb(int v, double j, double k) const;  // k1 - k2 is considered as a third parameter, k1 = ka, k2 = kc, 
	h2co_diagram(const std::string& path, molecule m, int& n_l, int verbosity = 1);
};


// The class contains the table with radiative coefficients;
// arr[i][j] is for i->j radiative transition; there is relation for rates: i->j *g_i = j->i *g_j
class einstein_coeff
{
protected:
	const energy_diagram *diagram;

public:
	int		nb_lev;
	double	**arr;

	einstein_coeff(const energy_diagram *di);
	virtual ~einstein_coeff();
};

// H2 molecule (Wolniewicz et al., Astroph. J. Suppl. Ser. 115, p. 293, 1998); (must be checked)
// only first 298 levels have one of the radiative transitions calculated by Wolniewicz et al. (1998);
class h2_einstein_coeff : public einstein_coeff
{
public:
	h2_einstein_coeff(const std::string &path, const energy_diagram *di, int verbosity =1);
};

// H2O molecule (HITRAN2012 database);
class h2o_einstein_coeff : public einstein_coeff
{
public:
	h2o_einstein_coeff(const std::string &path, const h2o_diagram *di, int verbosity =1);
};

// CH3OH molecule (Mekhtiev et al., Journal of Molecular Spectroscopy 194, 171-178, 1999);
class ch3oh_einstein_coeff : public einstein_coeff
{
public:
	ch3oh_einstein_coeff(const std::string &path, const energy_diagram *di, int verbosity =1);
};

// Radiative transitions of specific ion;
class ion_einstein_coeff : public einstein_coeff
{
public:
	ion_einstein_coeff(const std::string &path, const energy_diagram *, int verbosity =1);
};

// CO molecule radiative transitions;
class co_einstein_coeff : public einstein_coeff
{
public:
	co_einstein_coeff(const std::string &path, const energy_diagram *, int verbosity =1);
};

// OH molecule radiative transitions;
class oh_einstein_coeff : public einstein_coeff
{
public:
	oh_einstein_coeff(const std::string &path, const energy_diagram *, int verbosity =1);
};

// OH molecule radiative transitions, including hyperfine splitting;
class oh_hf_einstein_coeff : public einstein_coeff
{
public:
    oh_hf_einstein_coeff(const std::string& path, const energy_diagram*, int verbosity = 1);
};

// NH3 molecule radiative transitions;
class nh3_einstein_coeff : public einstein_coeff
{
public:
	nh3_einstein_coeff(const std::string &path, const energy_diagram *, int verbosity =1);
};

// H2CO molecule radiative transitions;
class h2co_einstein_coeff : public einstein_coeff
{
public:
	h2co_einstein_coeff(const std::string& path, const energy_diagram*, int verbosity = 1);
};


// The class containing the parameters of specific transition;
class transition
{
public:
	double energy, freq;
	energy_level low_lev, up_lev;

// The relation operators are needed to sort transition vector by energy:
	bool operator == (transition obj) const {return (low_lev == obj.low_lev && up_lev == obj.up_lev);}
	bool operator != (transition obj) const {return !(*this == obj);}
	bool operator < (transition obj) const {return (energy < obj.energy && !(*this == obj));}	
	bool operator > (transition obj) const {return (energy > obj.energy && !(*this == obj));}
	
	transition(energy_level lowl, energy_level upl);
};
