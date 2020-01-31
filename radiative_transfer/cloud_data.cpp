#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include "cloud_data.h"

using namespace std;
#define MAX_TEXT_LINE_WIDTH 6600 // must be long enough

cloud_layer::cloud_layer() :
	zl(0.), zu(0.), dz(0.), zm(0.), temp_n(0.), temp_el(0.), vel_n(0.), velg_n(0.), tot_h_conc(0.), he_conc(0.), h_conc(0.),
	oh2_conc(0.), ph2_conc(0.), el_conc(0.), mol_conc(0.), h2_opr(0.), vel_turb(0.), av_temp_d(0.)
{;}

cloud_layer::cloud_layer(const cloud_layer & obj) :
	zl(obj.zl), zu(obj.zu), dz(obj.dz), zm(obj.zm), 
	temp_n(obj.temp_n), temp_el(obj.temp_el), vel_n(obj.vel_n), velg_n(obj.velg_n), tot_h_conc(obj.tot_h_conc), he_conc(obj.he_conc), h_conc(obj.h_conc), 
	oh2_conc(obj.oh2_conc), ph2_conc(obj.ph2_conc), el_conc(obj.el_conc), mol_conc(obj.mol_conc), h2_opr(obj.h2_opr), vel_turb(obj.vel_turb), 
	av_temp_d(obj.av_temp_d)
{
	for (int i = 0; i < (int)obj.dust_grain_temp.size(); i++) {
		dust_grain_temp.push_back( obj.dust_grain_temp[i] );
	}

	for (int i = 0; i < (int)obj.dust_grain_conc.size(); i++) {
		dust_grain_conc.push_back( obj.dust_grain_conc[i] );
	}
}

cloud_layer::~cloud_layer() {
	dust_grain_temp.clear();
	dust_grain_conc.clear();
}

cloud_layer& cloud_layer::operator=(const cloud_layer& obj)
{
	zl = obj.zl; 
	zu = obj.zu;
	dz = obj.dz; 
	zm = obj.zm;
	
	temp_n = obj.temp_n;
	temp_el = obj.temp_el;
	vel_n = obj.vel_n;
	velg_n = obj.velg_n;
	tot_h_conc = obj.tot_h_conc; 
	he_conc = obj.he_conc; 
	h_conc = obj.h_conc;
	oh2_conc = obj.oh2_conc;
	ph2_conc = obj.ph2_conc;
	el_conc = obj.el_conc;
	mol_conc = obj.mol_conc; 
	h2_opr = obj.h2_opr; 
	vel_turb = obj.vel_turb;
	av_temp_d = obj.av_temp_d;

	dust_grain_temp.clear();
	for (int i = 0; i < (int)obj.dust_grain_temp.size(); i++) {
		dust_grain_temp.push_back(obj.dust_grain_temp[i]);
	}

	dust_grain_conc.clear();
	for (int i = 0; i < (int)obj.dust_grain_conc.size(); i++) {
		dust_grain_conc.push_back(obj.dust_grain_conc[i]);
	}
	return *this;
}


cloud_data::cloud_data() : nb_lay(0), dust(0)
{;}

void cloud_data::add_layer(cloud_layer & clayer) {
	lay_array.push_back(clayer);
	nb_lay++;
}

void cloud_data::remove_layer(int i)
{
	if (i > 0 && i < nb_lay) {
		lay_array.erase(lay_array.begin() + i);
		nb_lay--;
	}
}

void cloud_data::delete_layers() {
	lay_array.clear();
}

void cloud_data::set_vel_turb(double vt)
{
	for (int i = 0; i < (int) lay_array.size(); i++) {
		lay_array[i].set_vel_turb(vt);
	}
}

double cloud_data::get_height() const {
	return (lay_array.back().zu - lay_array.front().zl);
}

void join_layers(cloud_data* cloud, int nb)
{
	int i, j, l;
	double a;
	double* x = new double[nb];
	for (i = 0; i < nb*(cloud->nb_lay/nb); i += nb) 
	{
		// calculating weight for averaging,
        // weight of each layer is proportional to the H total column density,
		a = 0.;
		for (j = 0; j < nb; j++) {
			x[j] = cloud->lay_array[i+j].dz * cloud->lay_array[i+j].tot_h_conc;
			a += x[j];
		}
		for (j = 0; j < nb; j++) {
			x[j] /= a;
		}

		cloud->lay_array[i].zu = cloud->lay_array[i + nb-1].zu;
		cloud->lay_array[i].dz = cloud->lay_array[i].zu - cloud->lay_array[i].zl;
		cloud->lay_array[i].zm = cloud->lay_array[i].zl + 0.5 * cloud->lay_array[i].dz;
		
		cloud->lay_array[i].temp_n *= x[0];
		cloud->lay_array[i].temp_el *= x[0];
		cloud->lay_array[i].av_temp_d *= x[0];
		cloud->lay_array[i].vel_n *= x[0];
		cloud->lay_array[i].velg_n *= x[0];
		cloud->lay_array[i].tot_h_conc *= x[0];
		cloud->lay_array[i].he_conc *= x[0];
		cloud->lay_array[i].h_conc *= x[0];
		cloud->lay_array[i].oh2_conc *= x[0];
		cloud->lay_array[i].ph2_conc *= x[0];
		cloud->lay_array[i].el_conc *= x[0];
		cloud->lay_array[i].mol_conc *= x[0];
		cloud->lay_array[i].h2_opr *= x[0];
		cloud->lay_array[i].vel_turb *= x[0];

		for (l = 0; l < (int)cloud->lay_array[i].dust_grain_temp.size(); l++) {
			cloud->lay_array[i].dust_grain_temp[l] *= x[0];
			cloud->lay_array[i].dust_grain_conc[l] *= x[0];
		}

		for (j = 1; j < nb; j++) {
			cloud->lay_array[i].temp_n += cloud->lay_array[i + j].temp_n * x[j];
			cloud->lay_array[i].temp_el += cloud->lay_array[i + j].temp_el * x[j];
			cloud->lay_array[i].av_temp_d += cloud->lay_array[i + j].av_temp_d *x[j];
			cloud->lay_array[i].vel_n += cloud->lay_array[i + j].vel_n * x[j];
			cloud->lay_array[i].velg_n += cloud->lay_array[i + j].velg_n * x[j];
			cloud->lay_array[i].tot_h_conc += cloud->lay_array[i + j].tot_h_conc * x[j];
			cloud->lay_array[i].he_conc += cloud->lay_array[i + j].he_conc * x[j];
			cloud->lay_array[i].h_conc += cloud->lay_array[i + j].h_conc * x[j];
			cloud->lay_array[i].oh2_conc += cloud->lay_array[i + j].oh2_conc * x[j];
			cloud->lay_array[i].ph2_conc += cloud->lay_array[i + j].ph2_conc * x[j];
			cloud->lay_array[i].el_conc += cloud->lay_array[i + j].el_conc * x[j];
			cloud->lay_array[i].mol_conc += cloud->lay_array[i + j].mol_conc * x[j];
			cloud->lay_array[i].h2_opr += cloud->lay_array[i + j].h2_opr * x[j];
			cloud->lay_array[i].vel_turb += cloud->lay_array[i + j].vel_turb * x[j];
			
			for (l = 0; l < (int)cloud->lay_array[i].dust_grain_temp.size(); l++) {
				cloud->lay_array[i].dust_grain_temp[l] += cloud->lay_array[i + j].dust_grain_temp[l] * x[j];
				cloud->lay_array[i].dust_grain_conc[l] += cloud->lay_array[i + j].dust_grain_conc[l] * x[j];
			}
		}	
	}
	l = cloud->nb_lay/nb;
	for (i = 0; i < l; i++) {
		for (j = 1; j < nb; j++) {
			cloud->remove_layer(i + 1);
		}
	}
	
	j = (int) cloud->lay_array.size();
	for (i = l; i < j; i++) {
		cloud->remove_layer(l);
	}
}

void cloud_data::save_data(std::string fname)
{
    int i;
    ofstream output;

    output.open(fname.c_str());
    output << scientific;
  
    output << left << setw(18) << "!depth(cm)" << setw(13) << "gas_temp(K)" << setw(13) << "el_temp" << setw(13) << "dust_temp(K)"
        << setw(13) << "gasvel(cm/s)" << setw(13) << "concHe(cm-3)" << setw(13) << "conc_pH2" << setw(13) << "conc_oH2"
        << setw(13) << "conc_H" << setw(13) << "conc_e" << setw(13) << "conc_H_nucl" << setw(13) << "conc_mol" << endl;

    for (i = 0; i < nb_lay; i++) {
        output.precision(10);
        output << left << setw(18) << lay_array[i].zm;

        output.precision(4);
        output << left << setw(13) << lay_array[i].temp_n
            << setw(13) << lay_array[i].temp_el
            << setw(13) << lay_array[i].av_temp_d
            << setw(13) << lay_array[i].vel_n
            << setw(13) << lay_array[i].he_conc
            << setw(13) << lay_array[i].ph2_conc
            << setw(13) << lay_array[i].oh2_conc
            << setw(13) << lay_array[i].h_conc
            << setw(13) << lay_array[i].el_conc
            << setw(13) << lay_array[i].tot_h_conc
            << setw(13) << lay_array[i].mol_conc << endl;
    }
    output.close();
}


bool set_physical_parameters(std::string data_path, cloud_data* cloud)
{
    char text_line[MAX_TEXT_LINE_WIDTH];
    int i, j;
    double a, z, nh, ah, ah2, ahe, opr, td, abund;

    string fname;
    ifstream input;
    stringstream ss;
    cloud_layer clayer;

    // all files have empty line at the end,
    fname = data_path + "sim_phys_param.txt";
    input.open(fname.c_str());
    if (!input) {
        cout << "Cannot open file " << fname << endl;
        return false;
    }
    else {
        while (!input.eof())
        {
            do // comment lines are read:
                input.getline(text_line, MAX_TEXT_LINE_WIDTH);
            while (text_line[0] == '!' || text_line[0] == '#');

            if (text_line[0] == '\0')
                break;

            ss.clear();
            ss.str(text_line);
            // there are two velocity gradients - first is average, second is instantaneous from MHD equations
            ss >> clayer.zl >> a >> clayer.temp_n >> a >> clayer.temp_el >> clayer.vel_n >> a >> clayer.tot_h_conc >> a >> clayer.el_conc >> a
                >> a >> clayer.velg_n; // MHD velocity gradient

            clayer.el_conc *= clayer.tot_h_conc;
            cloud->add_layer(clayer);
        }
    }
    input.close();

    fname = data_path + "sim_data_h2_chemistry.txt";
    input.open(fname.c_str());
    if (!input) {
        cout << "Cannot open file " << fname << endl;
        return false;
    }
    else {
        i = 0;
        while (!input.eof())
        {
            do // comment lines are read:
                input.getline(text_line, MAX_TEXT_LINE_WIDTH);
            while (text_line[0] == '!' || text_line[0] == '#');

            if (text_line[0] == '\0')
                break;

            ss.clear();
            ss.str(text_line);

            ss >> z >> a;
            cloud->lay_array[i].h2_opr = a;
            i++;
        }
    }
    input.close();

    fname = data_path + "sim_specimen_abund.txt";
    input.open(fname.c_str());
    if (!input) {
        cout << "Cannot open file " << fname << endl;
        return false;
    }
    else {
        i = 0;
        while (!input.eof())
        {
            do // comment lines are read:
                input.getline(text_line, MAX_TEXT_LINE_WIDTH);
            while (text_line[0] == '!' || text_line[0] == '#');

            if (text_line[0] == '\0')
                break;

            ss.clear();
            ss.str(text_line);

            ss >> z >> ah >> ah2 >> ahe;

            nh = cloud->lay_array[i].tot_h_conc;
            opr = cloud->lay_array[i].h2_opr;

            cloud->lay_array[i].he_conc = ahe * nh;
            cloud->lay_array[i].ph2_conc = ah2 * nh / (1. + opr);
            cloud->lay_array[i].oh2_conc = ah2 * nh * opr / (1. + opr);
            cloud->lay_array[i].h_conc = ah * nh;
            i++;
        }
    }
    input.close();

    fname = data_path + "sim_dust_data.txt";
    input.open(fname.c_str());
    if (!input) {
        cout << "Cannot open file " << fname << endl;
        return false;
    }
    else {
        i = 0;
        while (!input.eof())
        {
            do // comment lines are read:
                input.getline(text_line, MAX_TEXT_LINE_WIDTH);
            while (text_line[0] == '!' || text_line[0] == '#');

            if (text_line[0] == '\0')
                break;

            ss.clear();
            ss.str(text_line);
            ss >> z;

            nh = cloud->lay_array[i].tot_h_conc;
            cloud->lay_array[i].dust_grain_temp.clear();
            cloud->lay_array[i].dust_grain_conc.clear();

            while (!ss.eof()) {
                ss >> td >> abund;
                for (j = 0; j < 14; j++) {
                    ss >> a;
                }
                cloud->lay_array[i].dust_grain_temp.push_back(td);
                cloud->lay_array[i].dust_grain_conc.push_back(abund * nh);
            }
            cloud->lay_array[i].av_temp_d = cloud->lay_array[i].dust_grain_temp.back();
            cloud->lay_array[i].dust_grain_temp.pop_back(); // removing average temperature and total abundance (charge)
            cloud->lay_array[i].dust_grain_conc.pop_back();
            i++;
        }
    }
    input.close();

    // calculating the average of physical parameters in the layer
    for (i = 0; i < cloud->nb_lay - 1; i++) {
        cloud->lay_array[i].zu = cloud->lay_array[i + 1].zl;
        cloud->lay_array[i].dz = cloud->lay_array[i].zu - cloud->lay_array[i].zl;
        cloud->lay_array[i].zm = cloud->lay_array[i].zl + 0.5 * cloud->lay_array[i].dz;

        cloud->lay_array[i].temp_n =     0.5*(cloud->lay_array[i].temp_n +     cloud->lay_array[i+1].temp_n);
        cloud->lay_array[i].temp_el =    0.5*(cloud->lay_array[i].temp_el +    cloud->lay_array[i+1].temp_el);
        cloud->lay_array[i].av_temp_d =  0.5*(cloud->lay_array[i].av_temp_d +  cloud->lay_array[i+1].av_temp_d);
        cloud->lay_array[i].vel_n =      0.5*(cloud->lay_array[i].vel_n +      cloud->lay_array[i+1].vel_n);
        cloud->lay_array[i].velg_n =     0.5*(cloud->lay_array[i].velg_n +     cloud->lay_array[i+1].velg_n);
        cloud->lay_array[i].tot_h_conc = 0.5*(cloud->lay_array[i].tot_h_conc + cloud->lay_array[i+1].tot_h_conc);
        cloud->lay_array[i].he_conc =    0.5*(cloud->lay_array[i].he_conc +    cloud->lay_array[i+1].he_conc);
        cloud->lay_array[i].h_conc =     0.5*(cloud->lay_array[i].h_conc +     cloud->lay_array[i+1].h_conc);
        cloud->lay_array[i].oh2_conc =   0.5*(cloud->lay_array[i].oh2_conc +   cloud->lay_array[i+1].oh2_conc);
        cloud->lay_array[i].ph2_conc =   0.5*(cloud->lay_array[i].ph2_conc +   cloud->lay_array[i+1].ph2_conc);
        cloud->lay_array[i].el_conc =    0.5*(cloud->lay_array[i].el_conc +    cloud->lay_array[i+1].el_conc);
        cloud->lay_array[i].mol_conc =   0.5*(cloud->lay_array[i].mol_conc +   cloud->lay_array[i+1].mol_conc);
        cloud->lay_array[i].h2_opr =     0.5*(cloud->lay_array[i].h2_opr +     cloud->lay_array[i+1].h2_opr);
        cloud->lay_array[i].vel_turb =   0.5*(cloud->lay_array[i].vel_turb +   cloud->lay_array[i+1].vel_turb);

        for (j = 0; j < (int)cloud->lay_array[i].dust_grain_temp.size(); j++) {
            cloud->lay_array[i].dust_grain_temp[j] = 0.5*(cloud->lay_array[i].dust_grain_temp[j] + cloud->lay_array[i+1].dust_grain_temp[j]);
            cloud->lay_array[i].dust_grain_conc[j] = 0.5*(cloud->lay_array[i].dust_grain_conc[j] + cloud->lay_array[i+1].dust_grain_conc[j]);
        }
    }
    cloud->remove_layer(cloud->nb_lay - 1);
    return true;
}

bool set_molecular_conc(std::string data_path, std::string mol_name, cloud_data* cloud, double f)
{
    char text_line[MAX_TEXT_LINE_WIDTH];
    int mol_nb, i, j, k;
    double a;

    string fname, str;
    vector<double> z_coord, abund, toth_conc;
    vector<double>::iterator it_low, it_up;

    ifstream input;
    stringstream ss;

    fname = data_path + "sim_specimen_abund.txt";
    input.open(fname.c_str());
    if (!input) {
        cout << "Cannot open file " << fname << endl;
        return false;
    }
    else {
        input.getline(text_line, MAX_TEXT_LINE_WIDTH);
        input.getline(text_line, MAX_TEXT_LINE_WIDTH);

        ss.clear();
        ss.str(text_line);
        ss >> str;

        mol_nb = 1;
        while (ss >> str) {
            if (str == mol_name)
                break;
            mol_nb++;
        }

        while (!input.eof())
        {
            input.getline(text_line, MAX_TEXT_LINE_WIDTH);
            if (text_line[0] == '\0')
                break;

            ss.clear();
            ss.str(text_line);

            ss >> a;
            z_coord.push_back(a);

            for (j = 0; j < mol_nb; j++) {
                ss >> a;
            }
            abund.push_back(a);
        }
    }
    input.close();

    fname = data_path + "sim_phys_param.txt";
    input.open(fname.c_str());
    if (!input) {
        cout << "Cannot open file " << fname << endl;
        return false;
    }
    else {
        while (!input.eof())
        {
            do // comment lines are read:
                input.getline(text_line, MAX_TEXT_LINE_WIDTH);
            while (text_line[0] == '!' || text_line[0] == '#');

            if (text_line[0] == '\0')
                break;

            ss.clear();
            ss.str(text_line);
            ss >> a >> a >> a >> a >> a >> a >> a >> a;
            toth_conc.push_back(a);
        }
    }
    input.close();

    for (i = 0; i < cloud->nb_lay; i++)
    {
        for (j = 0; j < (int)z_coord.size() - 1 && z_coord[j] < cloud->lay_array[i].zl; j++) { ; } // 		
        for (k = j; k < (int)z_coord.size() - 1 && z_coord[k] < cloud->lay_array[i].zu; k++) { ; } //

        if (j > 0 && z_coord[j] > cloud->lay_array[i].zl)
            cloud->lay_array[i].mol_conc += (z_coord[j] - cloud->lay_array[i].zl) * abund[j - 1] * toth_conc[j - 1];

        for (; j < k; j++) {
            cloud->lay_array[i].mol_conc += (z_coord[j + 1] - z_coord[j]) * abund[j] * toth_conc[j];
        }

        if (k > 0 && z_coord[k] > cloud->lay_array[i].zu)
            cloud->lay_array[i].mol_conc -= (z_coord[k] - cloud->lay_array[i].zu) * abund[k - 1] * toth_conc[k - 1];

        cloud->lay_array[i].mol_conc *= f / (cloud->lay_array[i].zu - cloud->lay_array[i].zl);
    }
    return true;
}
