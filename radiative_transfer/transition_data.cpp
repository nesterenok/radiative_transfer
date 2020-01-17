
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "integration.h"
#include "interpolation.h"
#include "transition_data.h"
#include "spectroscopy.h"

#define SOURCE_NAME "transition_data.cpp"
using namespace std;

transition_data::transition_data(int nb_lay, const energy_level &low, const energy_level &up)
: inv(0.), gain(0.), lum(0.), tau(0.), exc_temp(0.), loss_rate(0.), tau_sat(0.)
{
	nb_cloud_lay = nb_lay;
	trans = new transition(low, up);
	
	inv_arr = new double [nb_cloud_lay];
	memset(inv_arr, 0, nb_cloud_lay *sizeof(double));

	gain_arr = new double [nb_cloud_lay];
	memset(gain_arr, 0, nb_cloud_lay *sizeof(double));

	lum_arr = new double [nb_cloud_lay];
	memset(lum_arr, 0, nb_cloud_lay *sizeof(double));
}

transition_data::transition_data(const transition_data &obj)
{
	trans = new transition(obj.trans->low_lev, obj.trans->up_lev);
	
	nb_cloud_lay = obj.nb_cloud_lay;
	inv = obj.inv;
	gain = obj.gain;
	lum = obj.lum;
	tau = obj.tau;
	tau_sat = obj.tau_sat;
	exc_temp = obj.exc_temp;
	loss_rate = obj.loss_rate;
	
	inv_arr = new double [nb_cloud_lay];
	memcpy(inv_arr, obj.inv_arr, nb_cloud_lay *sizeof(double));

	gain_arr = new double [nb_cloud_lay];
	memcpy(gain_arr, obj.gain_arr, nb_cloud_lay *sizeof(double));

	lum_arr = new double [nb_cloud_lay];
	memcpy(lum_arr, obj.lum_arr, nb_cloud_lay *sizeof(double));
}

transition_data& transition_data::operator = (const transition_data &obj)
{
	delete trans;
	trans = new transition(obj.trans->low_lev, obj.trans->up_lev);
	
	nb_cloud_lay = obj.nb_cloud_lay;
	inv = obj.inv;
	gain = obj.gain;
	lum = obj.lum;
	tau = obj.tau;
	tau_sat = obj.tau_sat;
	exc_temp = obj.exc_temp;
	loss_rate = obj.loss_rate;

	delete [] inv_arr;
	inv_arr = new double [nb_cloud_lay];
	memcpy(inv_arr, obj.inv_arr, nb_cloud_lay *sizeof(double));

	delete [] gain_arr;
	gain_arr = new double [nb_cloud_lay];
	memcpy(gain_arr, obj.gain_arr, nb_cloud_lay *sizeof(double));

	delete [] lum_arr;
	lum_arr = new double [nb_cloud_lay];
	memcpy(lum_arr, obj.lum_arr, nb_cloud_lay *sizeof(double));

	return *this;
}

transition_data::~transition_data()
{
	delete trans;
	delete [] inv_arr;
	delete [] gain_arr;
	delete [] lum_arr;
}


transition_data_container::transition_data_container(const cloud_data *cl, const energy_diagram *di, const einstein_coeff *ec) 
{
	int i;
	saturated_maser_func func;

	cloud = cl;
	diagram = di;
	einst_coeff = ec;

	nb_mol_lev = diagram->nb_lev;
	nb_cloud_lay = cloud->nb_lay;
	fixed_nb = nb_cloud_lay / 2;

	t_nb = 35;
	tau_arr = new double [t_nb];
	sm_arr = new double [t_nb];

	for (i = 0; i < t_nb; i++)
	{
		tau_arr[i] = i;
		func.tau = tau_arr[i];
		sm_arr[i] = qromb<saturated_maser_func>(func, -5., 5., 1.e-5);
	}
}

transition_data_container::~transition_data_container()
{
	delete [] tau_arr;
	delete [] sm_arr;
}

const transition_data* transition_data_container::get(const transition *t) const
{
	list<transition_data>::const_iterator it;

	it = data.begin();
	while (it != data.end()) {
		if (it->trans->low_lev.nb == t->low_lev.nb && it->trans->up_lev.nb == t->up_lev.nb)
			return &(*it);
		it++;
	}
	return 0;
}

void transition_data_container::calc_inv(transition_data &trans_data, double *level_pop)
{
	double low_pop, up_pop;
	trans_data.inv = 0.;

	for (int lay = 0; lay < nb_cloud_lay; lay++)
	{
		low_pop = level_pop[lay*nb_mol_lev + trans_data.trans->low_lev.nb];
		up_pop = level_pop[lay*nb_mol_lev + trans_data.trans->up_lev.nb];
		
		trans_data.inv_arr[lay] = up_pop/trans_data.trans->up_lev.g - low_pop/trans_data.trans->low_lev.g;
		trans_data.inv += trans_data.inv_arr[lay] *cloud->lay_array[lay].dz;
	}
	trans_data.inv /= cloud->get_height();
}

void transition_data_container::calc_exc_temp(transition_data& trans_data, double* level_pop, int lay_nb)
{
	double low_pop, up_pop;
	low_pop = level_pop[lay_nb * nb_mol_lev + trans_data.trans->low_lev.nb];
	up_pop = level_pop[lay_nb * nb_mol_lev + trans_data.trans->up_lev.nb];
	
	trans_data.exc_temp = CM_INVERSE_TO_KELVINS *trans_data.trans->energy
		/log((low_pop *trans_data.trans->up_lev.g)/(up_pop *trans_data.trans->low_lev.g));
}

// The population inversion must be calculated before the function call;
void transition_data_container::calc_gain(transition_data &trans_data, double *level_pop)
{
	bool h2o_22GHz_case = false;
	int lay;
	double line_gain, d_abs, energy, energy_th, vel;
	
	energy = trans_data.trans->energy;
	energy_th = energy *energy *energy;
	
	// special case: H2O 6_16 -> 5_23 transition, the hyperfine splitting of 22.2 GHz line is taken into account;
	// check water molecule name,
    if (diagram->mol.name == "oH2O" && diagram->mol.isotop == 1 && trans_data.trans->up_lev.v == 0 
		&& trans_data.trans->up_lev.j == 6 && trans_data.trans->up_lev.k1 == 1 && trans_data.trans->up_lev.k2 == 6 
		&& trans_data.trans->low_lev.j == 5 && trans_data.trans->low_lev.k1 == 2 && trans_data.trans->low_lev.k2 == 3) 
		h2o_22GHz_case = true;

	for (lay = 0; lay < nb_cloud_lay; lay++) {
		if (h2o_22GHz_case) {	
			vel = sqrt(pow(sqrt(2.*BOLTZMANN_CONSTANT*cloud->lay_array[lay].temp_n/diagram->mol.mass) + 5.e+4, 2.)
					+ cloud->lay_array[lay].vel_turb *cloud->lay_array[lay].vel_turb);
		}
		else {
			vel = pow(2.*BOLTZMANN_CONSTANT*cloud->lay_array[lay].temp_n /diagram->mol.mass 
				+ cloud->lay_array[lay].vel_turb *cloud->lay_array[lay].vel_turb, 0.5);
		}

		// gain in the maser line at the line center, in [cm-1],
		line_gain = trans_data.inv_arr[lay] *trans_data.trans->up_lev.g 
			*einst_coeff->arr[trans_data.trans->up_lev.nb][trans_data.trans->low_lev.nb] *ONEDIVBY_SQRT_PI 
			*cloud->lay_array[lay].mol_conc /(energy_th *EIGHT_PI *vel);
		
		d_abs = cloud->dust->absorption(energy, cloud->lay_array[lay].dust_grain_conc);
		trans_data.gain_arr[lay] = line_gain - d_abs;
	}
	
	trans_data.gain = trans_data.tau = 0.;
	for (lay = 0; lay < nb_cloud_lay; lay++) 
	{
		trans_data.gain += trans_data.gain_arr[lay] *cloud->lay_array[lay].dz;
		// only layers with positive gain are taken into account in calculations of optical depth;
		if (trans_data.gain_arr[lay] > 0.) 
			trans_data.tau += trans_data.gain_arr[lay] *cloud->lay_array[lay].dz;
	}
	trans_data.gain /= cloud->get_height();
}

void transition_data_container::find(double *level_pop, double rel_error)
{
	bool is_inverted;
	int i, j, lay;
	transition_data *tr_data = 0;
	
	// clearing the old data on transitions:
	data.clear();

	for (i = 1; i < nb_mol_lev; i++) {
		for (j = 0; j < i; j++) {
			if (einst_coeff->arr[i][j] != 0)
			{
				tr_data = new transition_data(nb_cloud_lay, diagram->lev_array[j], diagram->lev_array[i]);
				calc_inv(*tr_data, level_pop);
				
				is_inverted = false;
				for (lay = 0; lay < nb_cloud_lay; lay++) {
					// the inversion must be higher than population error;
					is_inverted = (tr_data->inv_arr[lay] *diagram->lev_array[i].g > rel_error *level_pop[i]);
					if (is_inverted) break;
				}
				
				if (is_inverted) {
					calc_gain(*tr_data, level_pop);
					// the condition on optical depth value: 
					if (tr_data->tau > 0.1) {
						calc_exc_temp(*tr_data, level_pop, fixed_nb);
						data.push_front(*tr_data);
					}
				}
				delete tr_data;
			}
		}
	}
}

void transition_data_container::add(list<transition>& trans_list, double* level_pop)
{
    bool is_found;
    transition_data* tr_data = 0;

    list<transition>::iterator it;
    list<transition_data>::iterator it_d;

    it = trans_list.begin();
    while (it != trans_list.end())
    {
        // The search of the transition in the list of inverted transitions:
        is_found = false;
        it_d = data.begin();
        while (it_d != data.end())
        {
            if (*it == *(it_d->trans)) {
                is_found = true;
                break;
            }
            it_d++;
        }
        if (!is_found)
        {
            tr_data = new transition_data(nb_cloud_lay, it->low_lev, it->up_lev);
            calc_inv(*tr_data, level_pop);
            calc_gain(*tr_data, level_pop);
            calc_exc_temp(*tr_data, level_pop, fixed_nb);

            data.push_front(*tr_data);
            delete tr_data;
        }
        it++;
    }
    data.sort();
}

void transition_data_container::calc_saturation_depth(double beaming_factor)
{
	int i;
	double intensity, source_function;
	list<transition_data>::iterator it;

	it = data.begin();
	while (it != data.end())
	{
		it->tau_sat = 0.;
		if (it->exc_temp < 0.) 
		{
			source_function = 1./(1. - exp(it->trans->energy *CM_INVERSE_TO_KELVINS/it->exc_temp));
			intensity = it->loss_rate
				/(einst_coeff->arr[it->trans->up_lev.nb][it->trans->low_lev.nb] *beaming_factor *source_function);
		
			locate_index(sm_arr, t_nb, intensity, i);
            if (i >= 0 && i < t_nb - 1)
                it->tau_sat = tau_arr[i] + (tau_arr[i + 1] - tau_arr[i]) * (intensity - sm_arr[i]) / (sm_arr[i + 1] - sm_arr[i]);
            else if (i >= t_nb - 1)
                it->tau_sat = tau_arr[t_nb - 1];
		}
		it++;
	}
}

void transition_data_container::save_full_data(const string & fname) const
{
	int i, p;
	double z;
	
	ofstream outfile;
	list<transition_data>::const_iterator it;
	
	outfile.open(fname.c_str());
	outfile.setf(ios::scientific);
	
	if (!outfile.is_open()) 
		cout << "Error in" << SOURCE_NAME << ": can't open file to write transition data;" << endl;
	else 
	{
		outfile.precision(4);
		outfile << "Molecule name: " << diagram->mol.name << endl;
		outfile << "Physical parameters at a given cloud nb: " << fixed_nb << endl;
		outfile << left 
			<< setw(20) << "velocity gradient "	<< cloud->lay_array[fixed_nb].velg_n << " cm/s/cm" << endl
			<< setw(20) << "cloud height "		<< cloud->get_height() << " cm" << endl
			<< setw(20) << "turb velocity "		<< cloud->lay_array[fixed_nb].vel_turb << " cm/s" << endl
			<< setw(20) << "gas temperature "	<< cloud->lay_array[fixed_nb].temp_n << " K" << endl
			<< setw(20) << "h concentration "	<< cloud->lay_array[fixed_nb].h_conc << " cm-3" << endl
			<< setw(20) << "h2 concentration "	<< cloud->lay_array[fixed_nb].ph2_conc + cloud->lay_array[fixed_nb].oh2_conc << " cm-3" << endl
			<< setw(20) << "molecule concentr "	<< cloud->lay_array[fixed_nb].mol_conc << " cm-3" << endl
			<< setw(20) << "nb of mol levels "	<< diagram->nb_lev << endl << endl;

		outfile << "Dust:" << endl;
		outfile << left
			<< setw(16) << "name " << cloud->dust->name << endl
			<< setw(16) << "dust mass per H" << cloud->dust->dust_mass_perH << endl
			<< setw(16) << "temperature (K) ";
		for (i = 0; i < (int)cloud->lay_array[fixed_nb].dust_grain_temp.size(); i++) {
			outfile << left << setw(16) << cloud->lay_array[fixed_nb].dust_grain_temp[i];
		}
		outfile << endl << endl;
		
		outfile << "External radiation field:" << endl;
		for (i = 0; i < (int) cloud->ext_rad_field.size(); i++)
		{
			outfile << left << setw(17) << cloud->ext_rad_field[i]->name << endl
				<< setw(17) << "temperature (K) " << cloud->ext_rad_field[i]->temperature << endl
				<< setw(17) << "dilution " << cloud->ext_rad_field[i]->get_dilution() << endl << endl;
		}
		if ( (int) cloud->ext_rad_field.size() == 0) 
            outfile << "No" << endl << endl;

		outfile << "The number of transitions	" << data.size() << endl 
			<< "v, j, k1, k2, sym inv, (F) (upper) -> v, j, k1, k2, sym inv, (F) (lower)" << endl << endl;

		it = data.begin();
		while (it != data.end())
		{
			outfile.unsetf(ios_base::floatfield);
            outfile << left << setw(5) << it->trans->up_lev.v << setw(5) << it->trans->up_lev.j
                << setw(5) << it->trans->up_lev.k1 << setw(5) << it->trans->up_lev.k2 << setw(5) << it->trans->up_lev.syminv;

            if (diagram->hyperfine_splitting)
                outfile << left << setw(5) << it->trans->up_lev.hf;
                
            outfile << left << setw(5) << " -> " << setw(5) << it->trans->low_lev.v << setw(5) << it->trans->low_lev.j << setw(5) << it->trans->low_lev.k1
                << setw(5) << it->trans->low_lev.k2 << setw(5) << it->trans->low_lev.syminv;
			
            if (diagram->hyperfine_splitting)
                outfile << left << setw(5) << it->trans->low_lev.hf;
            
            outfile << left << setw(5) << setw(15) << setprecision(8) << it->trans->freq << endl << endl;
			
			outfile.setf(ios::scientific);
			outfile.precision(3);
			outfile << left 
				<< setw(15) << "mean inv " << it->inv << endl
				<< setw(15) << "mean gain " << it->gain << " cm-1" << endl
				<< setw(15) << "optical depth " << it->tau << endl
				<< setw(15) << "mean lum " << it->lum << " ph/cm3/s" << endl
				<< "At fixed layer:" << endl
				<< setw(15) << "exc temp " << it->exc_temp << " K" << endl
				<< setw(15) << "loss rate " << it->loss_rate << " s-1" << endl 
				<< setw(15) << "satur tau " << it->tau_sat << endl << endl;
			
			outfile << left << setw(17) << "relative depth " 
				<< setw(13) << "inversion "
				<< setw(13) << "gain "
				<< setw(13) << "luminosity " << endl;
			
			for (i = 0; i < nb_cloud_lay; i++)
			{
				if (i == 0) z = 0;
				else if (i == nb_cloud_lay-1) z = cloud->lay_array[i].zu;
				else z = cloud->lay_array[i].zm;
				
				// setting the precision
				p = (int) fabs(log10(1.-z/cloud->get_height())) + 2;
				if (p < 3) p = 3;
				if (p > 9) p = 9;

				outfile << left << setprecision(p) << setw(17) << z/cloud->get_height();
				outfile << left << setprecision(3) 
					<< setw(13) << it->inv_arr[i] 
					<< setw(13) << it->gain_arr[i] 
					<< setw(13) << it->lum_arr[i] << endl;
			}
			outfile << endl;
			it++;
		}
		outfile.close();
	}
}

void transition_data_container::save_transition(const std::string& fname, const transition &trans) const
{
    int i, p;
    double z;
    ofstream outfile;
  
    const transition_data* tr_data = get(&trans);

    outfile.open(fname.c_str());
    outfile.setf(ios::scientific);

    if (!outfile.is_open())
        cout << "Error in" << SOURCE_NAME << ": can't open file to write transition data;" << endl;
    else
    {
        outfile << left << setw(17) << "!depth(cm)" << setw(14) << "inv(cm-3)" << setw(14) << "gain(cm-1)" << setw(14) << "lum(ph/cm3/s)" << endl;
        for (i = 0; i < nb_cloud_lay; i++) 
        {
            if (i == 0) z = 0;
            else if (i == nb_cloud_lay - 1) z = cloud->lay_array[i].zu;
            else z = cloud->lay_array[i].zm;
            
            outfile << left << setprecision(6) << setw(17) << z;
            outfile << left << setprecision(4)
                << setw(14) << tr_data->inv_arr[i]
                << setw(14) << tr_data->gain_arr[i]
                << setw(14) << tr_data->lum_arr[i] << endl;
        }
    }
    outfile.close();
}

void transition_data_container::save_short_data(int process_nb, const string & fname) const
{
	double param1;
	ofstream outfile;
	list<transition_data>::const_iterator it;
	
	outfile.open(fname.c_str(), ios_base::out | ios_base::app);
	outfile.setf(ios::scientific);
	
	if (!outfile.is_open()) 
        cout << "Error in" << SOURCE_NAME << ": can't open file to write transition data;" << endl;
	else 
	{
		// Note: only the parameter of the first radiation field is included;
		if (!cloud->ext_rad_field.empty()) 
			param1 = cloud->ext_rad_field[0]->get_dilution();
		else param1 = 0.;

		outfile.precision(4);
		outfile << left << setw(5) << process_nb 
			<< setw(15) << cloud->get_height() 
			<< setw(15) << cloud->lay_array[fixed_nb].velg_n
			<< setw(15) << cloud->lay_array[fixed_nb].temp_n
			<< setw(15) << cloud->lay_array[fixed_nb].tot_h_conc
			<< setw(15) << cloud->lay_array[fixed_nb].mol_conc
			<< setw(15) << cloud->dust->dust_mass_perH
			<< setw(15) << cloud->lay_array[fixed_nb].dust_grain_temp[0]
			<< setw(15) << param1
			<< setw(15) << data.size() << endl;

		it = data.begin();
		while (it != data.end())
		{
			outfile << left 
				<< setw(5) << it->trans->up_lev.v << setw(5) << it->trans->up_lev.j << setw(5) << it->trans->up_lev.k1 
				<< setw(5) << it->trans->up_lev.k2 << setw(5) << it->trans->low_lev.v << setw(5) << it->trans->low_lev.j 
				<< setw(5) << it->trans->low_lev.k1 << setw(5) << it->trans->low_lev.k2 << setw(16) << setprecision(8) << it->trans->freq 
				<< setw(15) << setprecision(4) << it->inv << setw(15) << it->gain << setw(15) << it->lum << setw(15) << it->tau 
				<< setw(13) << it->exc_temp << setw(13) << setprecision(2) << it->tau_sat << endl; 	
			it++;
		}
		outfile.close();
	}
}
