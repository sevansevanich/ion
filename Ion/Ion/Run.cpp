#include "Particle.h"
#include <random> // for GSCH
#include <ctime> // fo GSCH
#include <fstream>

double C_res(int, double);
double B_res(double);
double n(int, double);
double k(double);
double E(double);
double E_I(double);
double E_norm(double, double);
double T_norm(double);

const double omega = 2 * M_PI*c / (0.8*pow(10, -6)); //circular frequency of laser beam 

//number of charge states be one upper then atomic state/charge

#pragma region input
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
} data_n;

struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 27;
	int ion_l[26] = { 0,0,2,2,2,2,2,2,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,0,0};
	int ion_n[26] = { 4, 4, 3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,1,1};
	double ion_pot[26] = { 7.9024678, 16.19920, 30.651, 54.91, 75.0, 98.985, 124.98, 151.060, 233.6, 262.10, 290.9, 330.8, 361.0, 392.2, 456.2, 489.312, 1262.7, 1357.8, 1460, 1575.6, 1687.0, 1798.4, 1950.4, 2045.759, 8828.1875, 9277.6814};
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

	data.number_of_charge_states = data_fe.number_of_charge_states;
	data.charge_step = data_fe.charge_step;
	data.charge = data_fe.charge;
	for (int i = 0;i < data_fe.number_of_charge_states;i++)
	{
		data.initial_pot(i, data_fe.ion_pot[i]);
		data.initial_l(i, data_fe.ion_l[i]);
		//data.initial_n(i, data_fe.ion_n[i]);
	}

	int z0 = 0, z1 = 0;

	z0 = data.charge / data.charge_step; // 
	z1 = z0 + data.number_of_charge_states - 1; //

	double *A = new double[z1];
	double *B = new double[z1];
	double *C = new double[z1];

	double fil = 0;
	Particle atom[kol][kol]; // target

	for (int i = 0;i < kol;i++) //target initialization for first ionization, n for the first electron
		for (int j = 0;j < kol;j++)
		{
			atom[i][j].charge = data.charge;
			atom[i][j].z1 = z1;
			atom[i][j].dt = 0.1*pow(10, -15);
		};

	double W[26] = {};
	int choose = 0, choose2 = 0;
	double fil0 = 0, value = 0, a0 = 0;

	cout << "Choose field's data:" << endl;
	cout << "1 - E" << endl;
	cout << "2 - I" << endl;
	cout << "3 - a" << endl;
	cin >> choose;
	cout << "Choose polarization:" << endl << "1 - linear" << endl << "2 - circular" << endl;
	cin >> choose2;
	cout << "Write value: ";
	cin >> value;
	value = pow(10, 22);

	switch (choose)
	{
	case 1:
		fil0 = value;
		break;
	case 2:
		if (choose2 == 1) {
			a0 = 8.6*pow(10, -10)*0.8*sqrt(value);// / sqrt(2); //intensy in W/cm^2 for circular polirazat
			fil0 = E(a0); //v/m
		}
		else
		{
			a0 = 8.6*pow(10, -10)*0.8*sqrt(value) / sqrt(2); //intensy in W/cm^2 for circular polirazat
			fil0 = E(a0); //v/m
		}
		break;
	case 3:
		fil0 = E(value);
		break;
	}

	//double r[21] = {};
	//for (int q = 0;q <= 20;q++)
	//	r[q] = abs(cos(-omega*q*pow(10, -15)));

	double ratio[kol][kol] = {};

	//ofstream fout;
	//fout.open("Ion_out_ion_E.txt");
	//double tau = 0.0;
	//while(tau<=40.0) //tau - time interval from the start of counting
	//{
	//	fil = fil0*exp(- pow(tau - 20, 2) / (2 * 40))*abs(sin(tau));
	//	fout << fil << " " << tau << " " << endl;
	//	tau += 0.01;
	//}
	//fout.close();

	double n_i[26] = {}, keld[26] = {};
	for (int i = 0; i < 26;i++)
	{
		n_i[i] = n(i + 1, k(data.ion_pot[i]));
		keld[i] = sqrt(2 * data.ion_pot[i] * m_e*1.6*pow(10, -19))*(2 * M_PI*c / (0.8*pow(10, -6))) / (el*fil); // ionization potential in Joule
	}

	for (int i = z0 + 1;i <= z1;i++)// calculation coefficients
	{
		int t = i - 1;
		C[t] = C_res(data.ion_l[t], n(i, k(data.ion_pot[t])));
		B[t] = B_res(data.ion_pot[t]);
		A[t] = w_a*pow(k(data.ion_pot[t]), 2)*(2 * data.ion_l[t] + 1)*pow(2 * B[t], (2 * n(i, k(data.ion_pot[t])) - 1))*C[t] / 2; //we devide by 2 because one 2 in degree 2*B^(2n-1)
	};

	ofstream fout;
	fout.open("Ion_out_ion_ratio.txt");
	int sum19=0, sum20=0, sum21 = 0, sum22 = 0, sum23 = 0, sum24 = 0, sum25 = 0;
	double tau = 0.1; //tau - time interval from the start of counting

	while(tau<= duration_imp) 
	{
		//fil = fil0*abs(sin(-omega*tsu*pow(10,-15))); // plane wave
		//fil = fil0*exp(-pow(tsu - 20, 2) / 40)); //gaus
		fil = fil0*exp(-pow(tau - 20, 2) / (2 * duration_imp))*abs(cos(omega*tau));

		/*ofstream fout;
		fout.open("Ion_out_ion.txt");

		double a = 0.1, W_E_a[30] = {};
		for (int i = 0; i < 30; i++)
		{
			W_E_a[i] = A[5] * pow(E(a), -(2 * n(6, k(data.ion_pot[5])) - 1))*exp(-2 * B[5] / (3 * E(a)));
			fout << W_E_a[i] << " " << a << " " << 1 - exp(-W_E_a[i] * 0.2668*pow(10, -14)) << endl;
			a += 0.1;
		}
		fout.close();*/

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
					ratio[i][j] = atom[i][j].charge / el;
					if ((ratio[i][j] <= 19.1) && (ratio[i][j] > 18.9)) { sum19++; };
					if ((ratio[i][j] <= 20.1) && (ratio[i][j] > 19.9)) { sum20++; };
					if ((ratio[i][j] <= 21.1) && (ratio[i][j] > 20.9)) { sum21++; };
					if ((ratio[i][j] <= 22.1) && (ratio[i][j] > 21.9)) { sum22++; };
					if ((ratio[i][j] <= 23.1) && (ratio[i][j] > 22.9)) { sum23++; };
					if ((ratio[i][j] <= 24.1) && (ratio[i][j] > 23.9)) { sum24++; };
					if ((ratio[i][j] <= 25.1) && (ratio[i][j] > 24.9)) { sum25++; };
				}
				atom[i][j].dt += 0.1*pow(10, -15); //increase time interval in code (fs)
			};
		fout<< sum19 << " "<< sum20 <<" " << sum21 << " " << sum22 << " " << sum23 << " " << sum24 << " " << sum25 << " " << tau << " " << endl;
		tau += 0.1; //increase time interval in loop (step)
		sum19 = 0; sum20 = 0; sum21 = 0; sum22 = 0; sum23 = 0; sum24 = 0; sum25 = 0;
	}

	fout.close();

	/*ofstream fout;
	fout.open("Ion_out_ion_ratio.txt");
	int sum21 = 0, sum22 = 0, sum23 = 0, sum24 = 0, sum25 = 0;
	for (int i = 0; i < kol; i++)
	{
		for (int j = 0; j < kol; j++)
		{
			ratio[i][j] = atom[i][j].charge / el;
			if ((ratio[i][j] <= 21.1) && (ratio[i][j] > 20.9)) { sum21++; };
			if ((ratio[i][j] <= 22.1) && (ratio[i][j] > 21.9)) { sum22++; };
			if ((ratio[i][j] <= 23.1) && (ratio[i][j] > 22.9)) { sum23++; };
			if ((ratio[i][j] <= 24.1) && (ratio[i][j] > 23.9)) { sum24++; };
			if ((ratio[i][j] <= 25.1) && (ratio[i][j] > 24.9)) { sum25++; };
		}
		
	}
	fout.close();*/

	delete[] A, B, C; //clearing memory from dinamic massiv
}


double C_res(int l_in, double n)
{
	if (n < 1) { n = 1; }; // not sure 
	double l = n - 1;
	//return pow(2 * M_PI*n, (-1))*pow((2 * M_E / n), 2 * n); // for l*<<n* 
	return pow(2, 2 * n) / (n*tgamma(n + l + 1)*tgamma(n - l)); //universal
	//return pow(2 * M_PI*n, (-1))*(;//universal without gamma function
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
double E(double a)
{
	//return (3.2*a / 0.8 )*pow(10, 12);
	return (a*m_e*c*c*2*M_PI)/(el*0.8*pow(10,-6)); //for wavelength 0.8 mkm
}

double E_I(double I)
{
	return sqrt(I*2*pow(10,4)/(c*epsilon0));
}

double E_norm(double E, double lym)
{
	return E*el*lym / (m_e*c*c*2*M_PI);
}

double T_norm(double lym)
{
	return lym / 2 * M_PI*c;
}