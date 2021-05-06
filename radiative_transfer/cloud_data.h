#pragma once
#include <vector>
#include <string>
#include "dust_model.h"

class cloud_data;
// the data are averaged over the group of nb layers, the last n < nb layers are removed,
void join_layers(cloud_data *cloud, int nb); 

// Reading the data provided by the simulation of shock waves, the full data path must be given,
// all old data is deleted, initializes all the data except mol_conc and vel_turb,
bool set_physical_parameters(std::string data_path, cloud_data* cloud);

// the specimen concentration is multiplied by factor f (in the case of spin isomers, e.g. o-H2 and p-H2),
// may be called after the joining of the cloud layers,
bool set_molecular_conc(std::string data_path, std::string mol_name, cloud_data* cloud, double f = 1.);

// initialization of level populations from file, provided by the simulation of shock waves,
// if given nb_lev != nb of levels in the file - error, lev_popul has dimension nb_lay * nb_lev, 
bool set_level_pop(std::string fname, cloud_data* cloud, int nb_lev, double *lev_popul);

// The cloud layer class
class cloud_layer
{
public:
	// layer coordinates,
	double zl, zu, dz, zm;
	// physical parameters of the gas,
	double temp_n, temp_el, av_temp_d, vel_n, velg_n, tot_h_conc, he_conc, h_conc, oh2_conc, ph2_conc, el_conc, mol_conc, h2_opr, vel_turb;
	
	// dust may contain several components having their own size, material, 
    // check the dust model used in simulations
	std::vector<double>  dust_grain_temp, dust_grain_conc; // must have the same size 
	
	void set_vel_turb(double vt) { vel_turb = vt; }
	void set_gas_temp(double tn) { temp_n = tn; }
	void set_electron_temp(double te) { temp_el = te; }
	void set_gas_velocity(double vn) { vel_n = vn; }
	void set_veln_grad(double vg) { velg_n = vg; }
	
	void set_toth_conc(double toth) { tot_h_conc = toth; }
	void set_h_conc(double c) { h_conc = c; }
	void set_he_conc(double c) { he_conc = c; }
	void set_oh2_conc(double c) { oh2_conc = c; }
	void set_ph2_conc(double c) { ph2_conc = c; }
	void set_el_conc(double c) { el_conc = c; }
	void set_mol_conc(double c) { mol_conc = c; }

	cloud_layer& operator = (const cloud_layer& obj);

	cloud_layer();
	cloud_layer(const cloud_layer&);
	~cloud_layer();
};

class cloud_data
{
public:
	int nb_lay;  // nb of cloud layers	
	const dust_model* dust;
	std::vector<const radiation_field*> ext_rad_field;
	std::vector<cloud_layer> lay_array;

	void add_radiation_field(const radiation_field* erf) { ext_rad_field.push_back(erf); }
	void set_dust_model(const dust_model* d) { dust = d; }
	void set_vel_turb(double);

    // returns the difference in height between last and first layers,
	double get_height() const;

	void add_layer(cloud_layer & layer);
	void remove_layer(int i); // i = 0 is the first layer
	void delete_layers();
    void save_data(std::string fname);

	cloud_data();
	virtual ~cloud_data() { delete_layers(); }
};
