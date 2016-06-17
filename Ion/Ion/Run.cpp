#include "Particle.h"
#include <random> // for GSCH
#include <ctime> // fo GSCH
#include <fstream>

double C_res(int, double);
double B_res(double);
double n(int, double);
double k(double);
double E(double, double);
double E_I(double);
double E_norm(double, double);
double T_norm(double);

//number of charge states be one upper then atomic state/charge

struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 2;
	int ion_l[1] = { 0 };
	int ion_n[1] = { 1 };
	int ion_num_of_el[1] = { 1 };
	double ion_pot[1] = { 13.598434005136 };
} data_h;

struct 
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 3;
	int ion_l[2] = { 0,0 };
	int ion_num_of_el[2] = { 2,1 };
	double ion_pot[2] = { 24.587387936,54.41776311 };
} data_he;

struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 7;
	int ion_l[6] = { 1,1,0,0,0,0 };
	int ion_n[6] = { 2,2,2,2,1,1 };
	int ion_num_of_el[6] = { 2,1,2,1,2,1 };
	double ion_pot[6] = { 11.26,24.38,47.88,64.49,392.09,489.99 };
} data_c;

struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 8;
	int ion_l[7] = { 1,1,1,0,0,0,0 };
	int ion_n[7] = { 2,2,2,2,2,1,1 };
	int ion_num_of_el[7] = { 3,2,1,2,1,2,1 };
	double ion_pot[7] = { 14.53413,29.60125,47.4453,77.4735,97.89013,552.06731,667.04609 };
	double density = 0.808;
	double atom_mass = 14.01;
} data_n;

struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 27;
	int ion_l[26] = { 0,0,2,2,2,2,2,2,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,0,0};
	int ion_n[26] = { 4, 4, 3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,1,1};
	double ion_pot[26] = { 7.9024678, 16.19920, 30.651, 54.91, 75.0, 98.985, 124.98, 151.060, 233.6, 262.10, 290.9, 330.8, 361.0, 392.2, 456.2, 489.312, 1262.7, 1357.8, 1460, 1575.6, 1687.0, 1798.4, 1950.4, 2045.759, 8828.1875, 9277.6814};
	double density = 7.874;
} data_fe;

struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 14;
	int ion_l[13] = { 1,0,0,1,1,1,1,1,1,0,0,0,0 };
	int ion_n[13] = { 3,2,2,2,2,2,2,2,2,2,2,1,1 };
	double ion_pot[13] = { 5.985768, 18.82855, 28.44764, 119.9924, 153.825, 190.49, 241.76, 284.64, 330.21, 398.65, 442.005, 2085.97963, 2304.14 };
} data_al;

struct data_in
{
	double charge;
	double charge_step;
	int number_of_charge_states;
	int *ion_l = new int[50];
	int *ion_n = new int[50];
	double *ion_pot = new double[50];
	double density;
	double atom_mass;

	//function for initialization of struct 
	void initial_pot(int i, double value)
	{
		ion_pot[i] = value;
	}
	void initial_l(int i, double value)
	{
		ion_l[i] = value;
	}
	void initial_n(int i, double value)
	{
		ion_n[i] = value;
	}
	void clear()
	{
		delete[] ion_l, ion_pot, ion_n;
	}
};

#pragma end_region

using namespace std;

int main()
{
	data_in data;
	ofstream fout;

	data.number_of_charge_states = data_n.number_of_charge_states;
	data.charge_step = data_n.charge_step;
	data.density = data_n.density;
	data.atom_mass = data_n.atom_mass;
	data.charge = data_n.charge;
	for (int i = 0;i < data_n.number_of_charge_states;i++)
	{
		data.initial_pot(i, data_n.ion_pot[i]);
		data.initial_l(i, data_n.ion_l[i]);
	}

	int z0 = 0, z1 = 0;

	z0 = data.charge / data.charge_step; // 
	z1 = z0 + data.number_of_charge_states - 1; //

	double *A = new double[z1];
	double *B = new double[z1];
	double *C = new double[z1];

	double lambda = pow(10, -6);
	double value = 0;

	cout << "Enter wavelength in micrometers for impulse: ";
	cin >> value;
	lambda *= value;

	double fil = 0;
	Particle atom[kol][kol]; // target
	double omega = 2 * M_PI*c / lambda; //circular frequency of laser beam 


	for (int i = 0;i < kol;i++) //target initialization for first ionization, n for the first electron
		for (int j = 0;j < kol;j++)
		{
			atom[i][j].charge = data.charge;
			atom[i][j].z1 = z1;
			atom[i][j].dt = 0.1*pow(10, -15);
		};

	double W[26] = {};
	int choose = 0, choose2 = 0;
	double fil0 = 0, a0 = 0, keld[26] = {};

	for (int i = z0 + 1;i <= z1;i++)// calculation coefficients
	{
		int t = i - 1;
		C[t] = C_res(data.ion_l[t], n(i, k(data.ion_pot[t])));
		B[t] = B_res(data.ion_pot[t]);
		A[t] = w_a*pow(k(data.ion_pot[t]), 2)*(2 * data.ion_l[t] + 1)*pow(2 * B[t], (2 * n(i, k(data.ion_pot[t])) - 1))*C[t] / 2; //we devide by 2 because one 2 in degree 2*B^(2n-1)
	};

	cout << "Choose what are you want: " << endl << " - calculation of depedence (1);" << endl << " - start to simulation of interaction of a laser pulse with a target (2);" << endl;
	cin >> choose; 

	if (choose == 1)
	{
		cout << "Choose charge of atom: ";
		cin >> choose2;

		fout.open("Ion_out_ion.txt");
		double n_i[26] = {};
		
		//Calculation of parameter Keldisha
		for (int i = 0; i < 26;i++)
		{
			n_i[i] = n(i + 1, k(data.ion_pot[i]));
			keld[i] = sqrt(2 * data.ion_pot[i] * m_e*1.6*pow(10, -19))*(omega) / (el*fil); // ionization potential in Joule
		}

		//Calculation of dependence of the ionization of the  nirmalize amplitude electric field
		double a = 0.001, W_E_a[100] = {};
		for (int i = 0; i < 100; i++)
		{
			W_E_a[i] = A[choose2] * pow(E(a, lambda), -(2 * n(choose2, k(data.ion_pot[choose2])) - 1))*exp(-2 * B[choose2] / (3 * E(a,lambda)));
			fout << W_E_a[i] << " " << a << " " << 1 - exp(-W_E_a[i] * 0.2668*pow(10, -14)) << endl;
			a += 0.0001;
		}
		fout.close();
	}
	else
	{
		cout << "Choose input field's data:" << endl;
		cout << "1 - E" << endl;
		cout << "2 - I" << endl;
		cout << "3 - a" << endl;
		cin >> choose;
		cout << "Choose polarization of a laser pulse:" << endl << "1 - linear" << endl << "2 - circular" << endl;
		cin >> choose2;
		cout << "Write value: ";
		cin >> value;
		//value = pow(10, 22);

		switch (choose)
		{
		case 1:
			fil0 = value;
			break;
		case 2:
			if (choose2 == 1) {
				a0 = 8.6*pow(10, -10)*0.8*sqrt(value);// / sqrt(2); //intensy in W/cm^2 for circular polirazat
				fil0 = E(a0, lambda); //v/m
			}
			else
			{
				a0 = 8.6*pow(10, -10)*0.8*sqrt(value) / sqrt(2); //intensy in W/cm^2 for circular polirazat
				fil0 = E(a0, lambda); //v/m
			}
			break;
		case 3:
			fil0 = E(value, lambda);
			break;
		}

		fout.open("Ion_out_ion_ratio_dif_N.txt");
		double n_e = 0.0;
		int *krat = new int[30];
		double tau = 0; //tau - time interval from the start of counting

		while (tau <= duration_imp)
		{
			//fil = fil0*abs(sin(-omega*tsu*pow(10,-15))); // plane wave
			//fil = fil0*exp(-pow(tsu - 20, 2) / 40)); //gaus
			//fil = fil0*exp(-pow(tau - 20, 2) * 3.45 / (2 * duration_imp)); //FWHM
			fil = fil0*exp(-pow(tau - 20, 2) * 3.45 / (2 * duration_imp))*abs(cos(omega*tau)); //FWHM with fill

			for (int i = 0;i < data.number_of_charge_states; i++) //clear charge states for new time step
				krat[i] = 0;

			for (int i = 0;i < kol;i++)
				for (int j = 0;j < kol;j++)
				{
					if (atom[i][j].z_i < z1)
					{
						for (int q = 0; q < z1;q++) // calculate ionization rate for z_i<i<z1. DON'T FORGET give current n
						{
							atom[i][j].W_cal(A[q], B[q], fil, n(q + 1, k(data.ion_pot[q])), keld[q]); // for direct current; give function n z, which equal q+1, but this equation does't >z1!!!!!!
							atom[i][j].W_cal_AC(A[q], B[q], fil, n(q + 1, k(data.ion_pot[q])), keld[q]); //for alternating current 
						}
						atom[i][j].P();

						if (atom[i][j].j != 0)
						{
							atom[i][j].charge += atom[i][j].j*data.charge_step; //Atom's charge will change if degree ionization doesn't equal zero (j!=0)
							atom[i][j].j = 0;
						};
						atom[i][j].clear();

						atom[i][j].dt += 0.1*pow(10, -15); //increase time interval in code (fs)

						for (int c = 0; c < data.number_of_charge_states; c ++) // calculate number of ions with charge=c after one time step
							if (c == atom[i][j].z_i) { krat[c]++; };
					}
				}
				
			for (int i = 0; i < data.number_of_charge_states; i++)
				fout << krat[i] << " "; //output to file

			fout << tau << endl;

			if ((tau <= 5) || (tau >= 27)) { tau += 1; } //increase time interval in loop (step)
			else {
				if ((tau >= 8) && (tau <= 16)) { tau += 0.1; }
				else { tau += 0.5; };
			}
			
		};			
		for (int i = 0; i < data.number_of_charge_states; i++) //calculate electorns density
			n_e += i*(data.density*N_av / data.atom_mass)*((double)krat[i]/((double)kol*(double)kol));
		fout << n_e << endl;
		delete[] krat; //clearing memory from dinamic massiv krat
	}
	
	fout.close();
	delete[] A, B, C; //clearing memory from dinamic massiv
}


double C_res(int l_in, double n)
{
	if (n < 1) { n = 1; }; // not sure 
	double l = n - 1;
	//return pow(2 * M_PI*n, (-1))*pow((2 * M_E / n), 2 * n); // for l*<<n*  for ADK
	return pow(2, 2 * n) / (n*tgamma(n + l + 1)*tgamma(n - l)); //universal
}
double B_res(double I)
{
	return (pow(k(I), 3)*E_a);
}
double n(int z, double k) // where z is charge atom after ionization
{
	return z / k;
}
double k(double I)
{
	return sqrt(I / I_h);
}
double E(double a, double lambda)
{
	return (a*m_e*c*c*2*M_PI)/(el*lambda); //for wavelength 0.8 mkm
}

double E_I(double I)
{
	return sqrt(I*2*pow(10,4)/(c*epsilon0));
}

double E_norm(double E, double lambda)
{
	return E*el*lambda / (m_e*c*c*2*M_PI);
}