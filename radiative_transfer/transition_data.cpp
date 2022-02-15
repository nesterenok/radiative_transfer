
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cstring>

#include "integration.h"
#include "interpolation.h"
#include "transition_data.h"
#include "spectroscopy.h"

#define SOURCE_NAME "transition_data.cpp"
using namespace std;

transition_data::transition_data(int nb_lay, const energy_level &low, const energy_level &up)
: inv(0.), gain(0.), lum(0.), tau_max(0.), tau_eff(0.), tau_sat(0.), temp_sat(0.), lay_nb_hg(0), 
delta_apect_ratio(0.25), nb_aspect_ratio(37), nb_freq(300)
{
	nb_cloud_lay = nb_lay;
	trans = new transition(low, up);
	
	inv_arr = new double [nb_cloud_lay];
	memset(inv_arr, 0, nb_cloud_lay *sizeof(double));

	gain_arr = new double [nb_cloud_lay];
	memset(gain_arr, 0, nb_cloud_lay *sizeof(double));

	lum_arr = new double [nb_cloud_lay];
	memset(lum_arr, 0, nb_cloud_lay *sizeof(double));

	emiss_coeff_arr = new double[nb_cloud_lay];
	memset(emiss_coeff_arr, 0, nb_cloud_lay * sizeof(double));

	pump_rate_arr = new double[nb_cloud_lay];
	memset(pump_rate_arr, 0, nb_cloud_lay * sizeof(double));

    pump_eff_arr = new double[nb_cloud_lay];
    memset(pump_eff_arr, 0, nb_cloud_lay * sizeof(double));

    loss_rate_arr = new double[nb_cloud_lay];
    memset(loss_rate_arr, 0, nb_cloud_lay * sizeof(double));

    exc_temp_arr = new double[nb_cloud_lay];
    memset(exc_temp_arr, 0, nb_cloud_lay * sizeof(double));

	tau_vs_aspect_ratio = new double[nb_aspect_ratio];
	memset(tau_vs_aspect_ratio, 0, nb_aspect_ratio * sizeof(double));

	tau_vs_frequency = new double[nb_freq];
	memset(tau_vs_frequency, 0, nb_freq * sizeof(double));
}

transition_data::transition_data(const transition_data &obj)
{
	trans = new transition(obj.trans->low_lev, obj.trans->up_lev);
	
	nb_aspect_ratio = obj.nb_aspect_ratio;
	nb_freq = obj.nb_freq;
	nb_cloud_lay = obj.nb_cloud_lay;
	lay_nb_hg = obj.lay_nb_hg;

	inv = obj.inv;
	gain = obj.gain;
	lum = obj.lum;
	tau_max = obj.tau_max;
	tau_eff = obj.tau_eff;
	tau_sat = obj.tau_sat;
	temp_sat = obj.temp_sat;
	delta_apect_ratio = obj.delta_apect_ratio;

	inv_arr = new double [nb_cloud_lay];
	memcpy(inv_arr, obj.inv_arr, nb_cloud_lay *sizeof(double));

	gain_arr = new double [nb_cloud_lay];
	memcpy(gain_arr, obj.gain_arr, nb_cloud_lay *sizeof(double));

	lum_arr = new double [nb_cloud_lay];
	memcpy(lum_arr, obj.lum_arr, nb_cloud_lay *sizeof(double));

	emiss_coeff_arr = new double[nb_cloud_lay];
	memcpy(emiss_coeff_arr, obj.emiss_coeff_arr, nb_cloud_lay * sizeof(double));

	pump_rate_arr = new double[nb_cloud_lay];
	memcpy(pump_rate_arr, obj.pump_rate_arr, nb_cloud_lay * sizeof(double));

    pump_eff_arr = new double[nb_cloud_lay];
    memcpy(pump_eff_arr, obj.pump_eff_arr, nb_cloud_lay * sizeof(double));

    loss_rate_arr = new double[nb_cloud_lay];
    memcpy(loss_rate_arr, obj.loss_rate_arr, nb_cloud_lay * sizeof(double));

    exc_temp_arr = new double[nb_cloud_lay];
    memcpy(exc_temp_arr, obj.exc_temp_arr, nb_cloud_lay * sizeof(double));

	tau_vs_aspect_ratio = new double[nb_aspect_ratio];
	memcpy(tau_vs_aspect_ratio, obj.tau_vs_aspect_ratio, nb_aspect_ratio * sizeof(double));

	tau_vs_frequency = new double[nb_freq];
	memcpy(tau_vs_frequency, obj.tau_vs_frequency, nb_freq * sizeof(double));
}

transition_data& transition_data::operator = (const transition_data &obj)
{
	delete trans;
	trans = new transition(obj.trans->low_lev, obj.trans->up_lev);
	
	nb_aspect_ratio = obj.nb_aspect_ratio;
	nb_freq = obj.nb_freq;
	nb_cloud_lay = obj.nb_cloud_lay;
	lay_nb_hg = obj.lay_nb_hg;

	inv = obj.inv;
	gain = obj.gain;
	lum = obj.lum;
	tau_max = obj.tau_max;
	tau_eff = obj.tau_eff;
	tau_sat = obj.tau_sat;
	temp_sat = obj.temp_sat;
	delta_apect_ratio = obj.delta_apect_ratio;
	
	delete [] inv_arr;
	inv_arr = new double [nb_cloud_lay];
	memcpy(inv_arr, obj.inv_arr, nb_cloud_lay *sizeof(double));

	delete [] gain_arr;
	gain_arr = new double [nb_cloud_lay];
	memcpy(gain_arr, obj.gain_arr, nb_cloud_lay *sizeof(double));

	delete [] lum_arr;
	lum_arr = new double [nb_cloud_lay];
	memcpy(lum_arr, obj.lum_arr, nb_cloud_lay *sizeof(double));

	delete[] emiss_coeff_arr;
	emiss_coeff_arr = new double[nb_cloud_lay];
	memcpy(emiss_coeff_arr, obj.emiss_coeff_arr, nb_cloud_lay * sizeof(double));

	delete[] pump_rate_arr;
	pump_rate_arr = new double[nb_cloud_lay];
	memcpy(pump_rate_arr, obj.pump_rate_arr, nb_cloud_lay * sizeof(double));

    delete[] pump_eff_arr;
    pump_eff_arr = new double[nb_cloud_lay];
    memcpy(pump_eff_arr, obj.pump_eff_arr, nb_cloud_lay * sizeof(double));

    delete[] loss_rate_arr;
    loss_rate_arr = new double[nb_cloud_lay];
    memcpy(loss_rate_arr, obj.loss_rate_arr, nb_cloud_lay * sizeof(double));

    delete[] exc_temp_arr;
    exc_temp_arr = new double[nb_cloud_lay];
    memcpy(exc_temp_arr, obj.exc_temp_arr, nb_cloud_lay * sizeof(double));

	delete[] tau_vs_aspect_ratio;
	tau_vs_aspect_ratio = new double[nb_aspect_ratio];
	memcpy(tau_vs_aspect_ratio, obj.tau_vs_aspect_ratio, nb_aspect_ratio * sizeof(double));

	delete[] tau_vs_frequency;
	tau_vs_frequency = new double[nb_freq];
	memcpy(tau_vs_frequency, obj.tau_vs_frequency, nb_freq * sizeof(double));

	return *this;
}

transition_data::~transition_data()
{
	delete trans;
	delete[] inv_arr;
	delete[] gain_arr;
    delete[] lum_arr;
	delete[] emiss_coeff_arr;
	delete[] pump_rate_arr;
    delete[] pump_eff_arr;
    delete[] loss_rate_arr;
    delete[] exc_temp_arr;
	delete[] tau_vs_aspect_ratio;
	delete[] tau_vs_frequency;
}


transition_data_container::transition_data_container(const cloud_data *cl, const energy_diagram *di, const einstein_coeff *ec) 
    : min_optical_depth(0.01)
{
	cloud = cl;
	diagram = di;
	einst_coeff = ec;

	nb_mol_lev = diagram->nb_lev;
	nb_cloud_lay = cloud->nb_lay;

	velocity_shift = 5.e+5;
}

transition_data_container::~transition_data_container()
{;}

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

	for (int lay = 0; lay < nb_cloud_lay; lay++) {
		low_pop = level_pop[lay*nb_mol_lev + trans_data.trans->low_lev.nb];
		up_pop = level_pop[lay*nb_mol_lev + trans_data.trans->up_lev.nb];
		
		trans_data.inv_arr[lay] = up_pop/trans_data.trans->up_lev.g - low_pop/trans_data.trans->low_lev.g;
		trans_data.inv += trans_data.inv_arr[lay] *cloud->lay_array[lay].dz;
	}
	trans_data.inv /= cloud->get_height();
}

void transition_data_container::calc_exc_temp(transition_data& trans_data, double* level_pop)
{
	int lay;
	double low_pop, up_pop;
    
	for (lay = 0; lay < nb_cloud_lay; lay++) {
        low_pop = level_pop[lay * nb_mol_lev + trans_data.trans->low_lev.nb];
        up_pop = level_pop[lay * nb_mol_lev + trans_data.trans->up_lev.nb];

        trans_data.exc_temp_arr[lay] = CM_INVERSE_TO_KELVINS * trans_data.trans->energy
            / log((low_pop * trans_data.trans->up_lev.g) / (up_pop * trans_data.trans->low_lev.g));
    }
}

// The population inversion must be calculated before the function call;
void transition_data_container::calc_gain(transition_data& trans_data)
{
	bool h2o_22GHz_case = false;
	int lay;
	double line_gain, d_abs, energy, energy_th, vel, g;

	energy = trans_data.trans->energy;
	energy_th = energy * energy * energy;

	// special case: H2O 6_16 -> 5_23 transition, the hyperfine splitting of 22.2 GHz line is taken into account;
	// check water molecule name,
	if (diagram->mol.name == "oH2O" && diagram->mol.isotop == 1 && trans_data.trans->up_lev.v == 0
		&& trans_data.trans->up_lev.j == 6 && trans_data.trans->up_lev.k1 == 1 && trans_data.trans->up_lev.k2 == 6
		&& trans_data.trans->low_lev.j == 5 && trans_data.trans->low_lev.k1 == 2 && trans_data.trans->low_lev.k2 == 3)
		h2o_22GHz_case = true;

	for (lay = 0; lay < nb_cloud_lay; lay++) {
		if (h2o_22GHz_case) {
			vel = sqrt(pow(sqrt(2. * BOLTZMANN_CONSTANT * cloud->lay_array[lay].temp_n / diagram->mol.mass) + 5.e+4, 2.)
				+ cloud->lay_array[lay].vel_turb * cloud->lay_array[lay].vel_turb);
		}
		else {
			vel = pow(2. * BOLTZMANN_CONSTANT * cloud->lay_array[lay].temp_n / diagram->mol.mass
				+ cloud->lay_array[lay].vel_turb * cloud->lay_array[lay].vel_turb, 0.5);
		}

		// gain in the maser line at the line center, in [cm-1],
		line_gain = trans_data.inv_arr[lay] * trans_data.trans->up_lev.g
			* einst_coeff->arr[trans_data.trans->up_lev.nb][trans_data.trans->low_lev.nb] * ONEDIVBY_SQRT_PI
			* cloud->lay_array[lay].mol_conc / (energy_th * EIGHT_PI * vel);

		d_abs = cloud->dust->absorption(energy, cloud->lay_array[lay].dust_grain_conc);
		trans_data.gain_arr[lay] = line_gain - d_abs;
	}

	trans_data.gain = trans_data.tau_eff = g = 0.;
	for (lay = 0; lay < nb_cloud_lay; lay++)
	{
		// gain may be negative here,
		trans_data.gain += trans_data.gain_arr[lay] * cloud->lay_array[lay].dz;

		// only layers with positive gain are taken into account in calculations of optical depth;
		if (trans_data.gain_arr[lay] > 0.)
			trans_data.tau_eff += trans_data.gain_arr[lay] * cloud->lay_array[lay].dz;

		if (trans_data.gain_arr[lay] > g) {
			g = trans_data.gain_arr[lay];
			trans_data.lay_nb_hg = lay;
		}
	}
	trans_data.gain /= cloud->get_height();

	if (g < 1.e-99)  // no inversion in the entire cloud, some very small value,
		trans_data.lay_nb_hg = 0;
}

void transition_data_container::calc_line_profile(transition_data & trans_data)
{
	int i, n, lay, nb_aspect_ratio, nb_freq;
	double da, x, vmin, vmax, dv, energy, energy_th, vel, profile, aspect_ratio;
	double* line_opacity, * vel_width, * dust_opacity;
	double** optical_depth;

	nb_freq = trans_data.nb_freq;
	nb_aspect_ratio = trans_data.nb_aspect_ratio;
	da = trans_data.delta_apect_ratio;

	energy = trans_data.trans->energy;
	energy_th = energy * energy * energy;

	// velocity range for spectra construction
	vmax = cloud->lay_array[0].vel_n + velocity_shift; // in cm/s
	vmin = cloud->lay_array[nb_cloud_lay - 1].vel_n - velocity_shift;
	dv = (vmax - vmin) / (nb_freq - 1.); 

	line_opacity = new double[nb_cloud_lay];
	memset(line_opacity, 0, nb_cloud_lay * sizeof(double));

	dust_opacity = new double[nb_cloud_lay];
	memset(dust_opacity, 0, nb_cloud_lay * sizeof(double));
	
	vel_width = new double[nb_cloud_lay];
	memset(vel_width, 0, nb_cloud_lay * sizeof(double));

	optical_depth = alloc_2d_array<double>(nb_freq, nb_aspect_ratio);
	memset(*optical_depth, 0, nb_freq * nb_aspect_ratio * sizeof(double));

	// H2O line 22 GHz need special treating (HF splitting expands the profile width)
	for (lay = 0; lay < nb_cloud_lay; lay++) {
		vel_width[lay] = pow(2. * BOLTZMANN_CONSTANT * cloud->lay_array[lay].temp_n / diagram->mol.mass
			+ cloud->lay_array[lay].vel_turb * cloud->lay_array[lay].vel_turb, 0.5);

		// inversion parameter is > 0 for population inversion, 
		line_opacity[lay] = trans_data.inv_arr[lay] * trans_data.trans->up_lev.g
			* einst_coeff->arr[trans_data.trans->up_lev.nb][trans_data.trans->low_lev.nb] * cloud->lay_array[lay].mol_conc *ONEDIVBY_SQRT_PI 
			/(energy_th * EIGHT_PI * vel_width[lay]);

		dust_opacity[lay] = cloud->dust->absorption(energy, cloud->lay_array[lay].dust_grain_conc);
	}

	for (i = 0; i < nb_aspect_ratio; i++) {
		aspect_ratio = 1. + da * i;  // = 1/mu
		
		for (n = 0, vel = vmin; n < nb_freq; vel += dv, n++) {
			optical_depth[n][i] = 0.;

			for (lay = 0; lay < nb_cloud_lay; lay++) {
				x = (vel - cloud->lay_array[lay].vel_n / aspect_ratio) / vel_width[lay];
				profile = exp(-x * x);

				// positive for population inversion
				if (line_opacity[lay] * profile - dust_opacity[lay] > 0.)
					optical_depth[n][i] += (line_opacity[lay] *profile - dust_opacity[lay]) 
						* cloud->lay_array[lay].dz *aspect_ratio;
			}
		}
	}

	// maximal absolute value of tau along the outflow
	for (i = 0; i < nb_aspect_ratio; i++) {
		x = 0.;
		for (n = 0; n < nb_freq; n++) {
			if (x < optical_depth[n][i])
				x = optical_depth[n][i];
		}
		trans_data.tau_vs_aspect_ratio[i] = x;
	}
	trans_data.tau_max = trans_data.tau_vs_aspect_ratio[0];
	
    // the dependence of optical depth on frequency for the direction along the gas flow,
	for (n = 0; n < nb_freq; n++) {
		trans_data.tau_vs_frequency[n] = optical_depth[n][0];
	}

	delete[] line_opacity; 
	delete[] vel_width; 
	delete[] dust_opacity;
	free_2d_array(optical_depth);
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
			if (einst_coeff->arr[i][j] != 0) {
				tr_data = new transition_data(nb_cloud_lay, diagram->lev_array[j], diagram->lev_array[i]);
				calc_inv(*tr_data, level_pop);
				
				is_inverted = false;
				for (lay = 0; lay < nb_cloud_lay; lay++) {
					// the inversion must be higher than population error;
					is_inverted = (tr_data->inv_arr[lay] *diagram->lev_array[i].g > rel_error *level_pop[i]);
					if (is_inverted) 
                        break;
				}
				
				if (is_inverted) {
					calc_gain(*tr_data);
					calc_line_profile(*tr_data);

					// the condition on optical depth value: 
					if (tr_data->tau_max >= min_optical_depth) {
						calc_exc_temp(*tr_data, level_pop);
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
    while (it != trans_list.end()) {
        // the search of the transition in the list of inverted transitions:
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
            calc_gain(*tr_data);
			calc_line_profile(*tr_data);
            calc_exc_temp(*tr_data, level_pop);
			
            data.push_front(*tr_data);
            delete tr_data;
        }
        it++;
    }
    data.sort();
}

void transition_data_container::calc_saturation_depth(double beaming_factor)
{
	int l;
	double source_function, cmb_intensity, exc_temp, intensity_0;
	list<transition_data>::iterator it;

	it = data.begin();
	while (it != data.end())
	{
		it->tau_sat = 0.;
		l = it->lay_nb_hg;
		exc_temp = it->exc_temp_arr[l];

        if (exc_temp < 0.) {
            source_function = 1. / (1. - exp(it->trans->energy * CM_INVERSE_TO_KELVINS / exc_temp));
			cmb_intensity = 1. / (exp(it->trans->energy * CM_INVERSE_TO_KELVINS / 2.73) - 1.);
			intensity_0 = source_function + cmb_intensity;
	
			// loss rate is an average for upper and lower levels,
			// without averaging over the line profile
			it->tau_sat = log(it->loss_rate_arr[l] 
				/ (einst_coeff->arr[it->trans->up_lev.nb][it->trans->low_lev.nb] * beaming_factor * intensity_0));
			
			it->temp_sat =  intensity_0 * exp(it->tau_sat) *PLANCK_CONSTANT * it->trans->freq / BOLTZMANN_CONSTANT;
        }       
		it++;
	}
}

void transition_data_container::save_full_data(const string & fname) const
{
	int i, p, fixed_nb;
	double z;
	
	ofstream outfile;
	list<transition_data>::const_iterator it;
	
	outfile.open(fname.c_str());
	outfile.setf(ios::scientific);
	
	if (!outfile.is_open()) 
		cout << "Error in " << SOURCE_NAME << ": can't open file to write transition data: " << endl << "	" << fname << endl;
	else 
	{ 
		// cloud parameters are saved for the cloud centre,
		fixed_nb = nb_cloud_lay / 2;

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
			
			// new definition of the fixed nb, 
			fixed_nb = it->lay_nb_hg;

			outfile.setf(ios::scientific);
			outfile.precision(3);
			outfile << left 
				<< setw(15) << "mean inv " << it->inv << endl
				<< setw(15) << "mean gain " << it->gain << " cm-1" << endl
				<< setw(15) << "optical depth (eff) " << it->tau_eff << endl
				<< setw(15) << "optical depth (line centre) " << it->tau_max << endl
				<< setw(15) << "mean lum " << it->lum << " ph/cm3/s" << endl
				<< "At layer with the highest gain:" << endl
				<< setw(15) << "inv " << it->inv_arr[fixed_nb] << endl
				<< setw(15) << "gain (cm-1)" << it->gain_arr[fixed_nb] << endl
				<< setw(15) << "lum (ph/cm3/s) " << it->lum_arr[fixed_nb] << endl
				<< setw(15) << "satur tau " << it->tau_sat << endl 
				<< setw(15) << "satur temp (K) " << it->temp_sat << endl << endl;
			
			outfile << left << setw(17) << "relative depth " 
				<< setw(12) << "inversion "
				<< setw(12) << "gain "
                << setw(12) << "exc temp"
				<< setw(12) << "luminosity "
				<< setw(12) << "emiss meas "
				<< setw(12) << "pump rate"
                << setw(12) << "pump effic "
                << setw(12) << "loss rate " << endl;
			
			for (i = 0; i < nb_cloud_lay; i++) {
				if (i == 0) z = 0;
				else if (i == nb_cloud_lay-1) z = cloud->lay_array[i].zu;
				else z = cloud->lay_array[i].zm;
				
				// setting the precision
				p = (int) fabs(log10(1. - z/cloud->get_height())) + 2;
				if (p < 3) p = 3;
				if (p > 9) p = 9;

				outfile << left << setprecision(p) << setw(17) << z/cloud->get_height();
				outfile << left << setprecision(3) 
					<< setw(12) << it->inv_arr[i] 
					<< setw(12) << it->gain_arr[i] 
                    << setw(12) << it->exc_temp_arr[i]
					<< setw(12) << it->lum_arr[i] 
					<< setw(12) << it->emiss_coeff_arr[i]
					<< setw(12) << it->pump_rate_arr[i]
                    << setw(12) << it->pump_eff_arr[i]
                    << setw(12) << it->loss_rate_arr[i] << endl;
			}
			outfile << endl;
			it++;
		}
		outfile.close();
	}
}

void transition_data_container::save_short_data(int process_nb, const string & fname) const
{
    int fixed_nb;
	double rf_param;
	ofstream outfile;
	list<transition_data>::const_iterator it;
	
	outfile.open(fname.c_str(), ios_base::out | ios_base::app);
	outfile.setf(ios::scientific);
	
	if (!outfile.is_open()) 
        cout << "Error in " << SOURCE_NAME << ": can't open file to write transition data: " << fname << endl;
	else {
		// cloud parameters are saved for the cloud centre,
        fixed_nb = nb_cloud_lay / 2;

		if (!cloud->ext_rad_field.empty()) // note: only the parameter of the first radiation field is included;
			rf_param = cloud->ext_rad_field[0]->get_dilution();
		else rf_param = 0.;
		
		outfile << "! Phys. parameters (at cloud centre): " 
			<< "height(cm), dv/dz(cm/s/cm), Tn(K), n_Htot(cm-3), conc_mol(cm-3), dust_mass_perH(g), Td(K), rad_field_parameter, nb" << endl
			<< "! Maser parameters (at layer with the highest gain): " 
			<< "v_u, j_u, k1_u, k2_u, sym_u, v_l, j_l, k1_l, k2_l, sym_l, freq(Hz), inversion(normalized_on_conc), gain(cm-1), luminosity(ph/cm3/s), opt_depth_max(along_the_gas_flow), T_exc(K), opt_depth_sat" << endl;

		outfile.setf(ios::scientific);
		outfile.precision(3);
		outfile << left << setw(5) << process_nb 
			<< setw(12) << cloud->get_height() 
			<< setw(12) << cloud->lay_array[fixed_nb].velg_n
			<< setw(12) << cloud->lay_array[fixed_nb].temp_n
			<< setw(12) << cloud->lay_array[fixed_nb].tot_h_conc
			<< setw(12) << cloud->lay_array[fixed_nb].mol_conc
			<< setw(12) << cloud->dust->dust_mass_perH
			<< setw(12) << cloud->lay_array[fixed_nb].dust_grain_temp[0]
			<< setw(12) << rf_param
			<< setw(12) << data.size() << endl;

		it = data.begin();
		while (it != data.end())
		{
			// new definition of the fixed nb,
			fixed_nb = it->lay_nb_hg;

			outfile.unsetf(ios_base::floatfield);
			outfile << left
				<< setw(5) << it->trans->up_lev.v << setw(5) << it->trans->up_lev.j << setw(5) << it->trans->up_lev.k1
				<< setw(5) << it->trans->up_lev.k2 << setw(5) << it->trans->up_lev.syminv 
				<< setw(5) << it->trans->low_lev.v << setw(5) << it->trans->low_lev.j << setw(5) << it->trans->low_lev.k1 
				<< setw(5) << it->trans->low_lev.k2 << setw(5) << it->trans->low_lev.syminv;
			
			outfile.setf(ios::scientific);
			outfile << setw(17) << setprecision(8) << it->trans->freq
				<< setw(12) << setprecision(3) << it->inv_arr[fixed_nb] << setw(12) << it->gain_arr[fixed_nb] << setw(12) << it->lum_arr[fixed_nb] << setw(12) << it->tau_max
				<< setw(12) << it->exc_temp_arr[fixed_nb] << setw(11) << setprecision(2) << it->tau_sat << endl;
			it++;
		}
		outfile.close();
	}
}

void transition_data_container::save_optical_depth_1(const std::string& fname) const
{
	int i, nb_aspect_ratio;
	double da;
	ofstream outfile;
	list<transition_data>::const_iterator it;

	outfile.open(fname.c_str());
	outfile.setf(ios::scientific);

	if (!outfile.is_open())
		cout << "Error in " << SOURCE_NAME << ": can't open file to write transition data: " << endl << "	" << fname << endl;
	else {
		outfile << left << "!" << endl << setw(14) << "!";
		it = data.begin();
		while (it != data.end()) {
			outfile << setw(14) << setprecision(6) << it->trans->freq;
			it++;
		}
		outfile << endl;

		nb_aspect_ratio = data.begin()->nb_aspect_ratio;
		da = data.begin()->delta_apect_ratio;
		
		for (i = 0; i < nb_aspect_ratio; i++) {
			outfile << left << setw(14) << setprecision(3) << 1. + da * i;
			it = data.begin();
			while (it != data.end()) {
				outfile << setw(14) << setprecision(3) << it->tau_vs_aspect_ratio[i];
				it++;
			}
			outfile << endl;
		}
	}
}

void transition_data_container::save_optical_depth_2(const std::string& fname) const
{
	int i;
	double dv, vmin, vmax;
	ofstream outfile;
	list<transition_data>::const_iterator it;

	outfile.open(fname.c_str());
	outfile.setf(ios::scientific);

	if (!outfile.is_open())
		cout << "Error in " << SOURCE_NAME << ": can't open file to write transition data: " << endl << "	" << fname << endl;
	else {
		outfile << left << "! gas is moving from the observer, reference frame of the shock front (see also Flower et al., MNRAS 409, 29, 2010)"
			<< endl << setw(14) << "! ";
		it = data.begin();
		while (it != data.end()) {
			outfile << setw(14) << setprecision(6) << it->trans->freq;
			it++;
		}
		outfile << endl;

		vmax = cloud->lay_array[0].vel_n + velocity_shift;
		vmin = cloud->lay_array[nb_cloud_lay - 1].vel_n - velocity_shift;
		dv = (vmax - vmin) / (data.begin()->nb_freq - 1.);

		for (i = 0; i < data.begin()->nb_freq; i++) {
			outfile << left << setw(14) << setprecision(4) << vmin + i*dv;  // cloud->lay_array[0].vel_n - ...
			it = data.begin();
			while (it != data.end()) {
				outfile << setw(14) << setprecision(3) << it->tau_vs_frequency[i];
				it++;
			}
			outfile << endl;
		}
	}
}

void transition_data_container::save_transition(const std::string& fname, const transition &trans) const
{
    int i;
    double z;
    ofstream outfile;
    const transition_data* tr_data = get(&trans);

    outfile.open(fname.c_str());
    outfile.setf(ios::scientific);

    if (!outfile.is_open())
        cout << "Error in " << SOURCE_NAME << ": can't open file to write transition data: " << fname << endl;
    else
    {
        outfile << left << setw(12) << "!depth(cm)" << setw(12) << "dz(cm)" << setw(12) << "inv(cm-3)" << setw(12) << "gain(cm-1)" 
            << setw(12) << "exctemp(K)" << setw(14) << "lum(ph/cm3/s)" << setw(12) << "emiss_meas" << setw(12) << "pump(cm3/s)" << setw(10) << "pumpeffic" 
			<< setw(14) << "lossrate(s-1)" << endl;
        
		for (i = 0; i < nb_cloud_lay; i++) {
			// for first and last layer - cloud boundary is set, for other layers - the middle point of the layer
            if (i == 0) z = 0;
            else if (i == nb_cloud_lay - 1) z = cloud->lay_array[i].zu;
            else z = cloud->lay_array[i].zm;
            
            // absolute depth and dz are saved (not normalized), 
            outfile << left << setprecision(3) << setw(12) << z << setw(12) << cloud->lay_array[i].dz
                << setw(12) << tr_data->inv_arr[i]
                << setw(12) << tr_data->gain_arr[i]
                << setw(12) << tr_data->exc_temp_arr[i]
                << setw(12) << tr_data->lum_arr[i] 
				<< setw(12) << tr_data->emiss_coeff_arr[i]
				<< setw(12) << tr_data->pump_rate_arr[i]
                << setw(12) << tr_data->pump_eff_arr[i]
                << setw(12) << tr_data->loss_rate_arr[i] << endl;
        }
    }
    outfile.close();
}
