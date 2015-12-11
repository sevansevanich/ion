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
	int ion_l[7] = { 1,1,1,0,0,0,0 };
	int ion_n[7] = { 2,2,2,2,2,1,1 };
	int ion_num_of_el[7] = { 3,2,1,2,1,2,1 };
	double ion_pot[7] = { 14.53413,29.60125,47.4453,77.4735,97.89013,552.06731,667.04609 };
} data_n;

struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 6;
	int ion_l[6] = { 1,1,0,0,0,0 };
	int ion_n[6] = { 2,2,2,2,1,1 };
	int ion_num_of_el[6] = { 2,1,2,1,2,1 };
	double ion_pot[6] = { 11.26,24.38,47.88,64.49,392.09,489.99 };
} data_c;

struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 26;
	int ion_l[26] = { 0,0,2,2,2,2,2,2,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,0,0};
	int ion_n[26] = { 4, 4, 3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,1,1};
	double ion_pot[26] = { 7.9024678, 16.19920, 30.651, 54.91, 75.0, 98.985, 124.98, 151.060, 233.6, 262.10, 290.9, 330.8, 361.0, 392.2, 456.2, 489.312, 1262.7, 1357.8, 1460, 1575.6, 1687.0, 1798.4, 1950.4, 2045.759, 8828.1875, 9277.6814};
} data_fe;

struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 13;
	int ion_l[13] = { 1,0,0,1,1,1,1,1,1,0,0,0,0 };
	int ion_n[13] = { 3,2,2,2,2,2,2,2,2,2,2,1,1 };
	double ion_pot[13] = { 5.985768, 18.82855, 28.44764, 119.9924, 153.825, 190.49, 241.76, 284.64, 330.21, 398.65, 442.005, 2085.97963, 2304.14 };
} data_al;

struct data_in
{
	double charge;
	double charge_step;
	int number_of_charge_states;
	int *ion_l = new int[25];
	int *ion_n = new int[25];
	int *ion_num_of_el = new int[25];
	double ion_pot[];
};

using namespace std;

int main()
{
	data_in data;

	int z0 = 0, z1 = 0;

	z0 = data_n.charge / data_n.charge_step; // 
	z1 = z0 + data_n.number_of_charge_states - 1; //

	double *A = new double[z1];
	double *B = new double[z1];
	double *C = new double[z1];

	double W[10] = {};
	int choose = 0, choose2 = 0;
	double fil = 0, value = 0, a0 = 0;

	cout << "Choose field's data:" << endl;
	cout << "1 - E" << endl;
	cout << "2 - I" << endl;
	cout << "3 - a" << endl;
	cin >> choose;
	cout << "Choose polarization:" << endl << "1 - linear" << endl << "2 - circular" << endl;
	cin >> choose2;
	cout << "Write value: ";
	cin >> value;
	//value = 9.23*pow(10, 13);

	switch (choose)
	{
	case 1:
		fil = value;
		break;
	case 2:
		if (choose2 == 1) {
			a0 = 8.6*pow(10, -10)*0.8*sqrt(value);// / sqrt(2); //intensy in W/cm^2 for circular polirazat
			fil = E(a0); //v/m
		}
		else
		{
			a0 = 8.6*pow(10, -10)*0.8*sqrt(value) / sqrt(2); //intensy in W/cm^2 for circular polirazat
			fil = E(a0); //v/m
		}
		break;
	case 3:
		fil = E(value);
		break;
	}

	for (int i = z0 + 1;i <= z1;i++)// calculation coefficients
	{
		int t = i - 1;
		C[t] = C_res(data_n.ion_l[t], n(i, k(data_n.ion_pot[t])));
		B[t] = B_res(data_n.ion_pot[t]);
		A[t] = w_a*pow(k(data_n.ion_pot[t]), 2)*(2 * data_n.ion_l[t] + 1)*pow(B[t] * 2, (2 * n(i, k(data_n.ion_pot[t])) - 1))*C[t] / 2;
		W[t] = A[t] * pow(E(2), -(2 * n(i, k(data_n.ion_pot[t])) - 1))*exp(-2 * B[t] / (3 * E(2)));
	};

	ofstream fout;
	fout.open("Ion_out_ion.txt");

	double a = 0.1, W_E_a[30] = {};
	for (int i = 0; i < 30; i++)
	{
		W_E_a[i] = A[5] * pow(E(a), -(2 * n(6, k(data_n.ion_pot[5])) - 1))*exp(-2 * B[5] / (3 * E(a)));
		fout << W_E_a[i] << " " << a << " " << 1 - exp(-W_E_a[i] * 0.2668*pow(10, -14)) << endl;
		a += 0.1;
	}
	fout.close();
	
	Particle atom[kol][kol]; // target

	for (int i = 0;i < kol;i++) //target initialization for first ionization, n for the first electron
		for (int j = 0;j < kol;j++)
		{
			atom[i][j].charge = data_n.charge;
			atom[i][j].z1 = z1;
			atom[i][j].n = n(1, k(data_n.ion_pot[0]));
			atom[i][j].dt = 0.2668*pow(10, -14);
		};
	
	for (int i = 0;i < kol;i++)
		for (int j = 0;j < kol;j++)
		{
			if (atom[i][j].z_i<z1)
			{
				for (int q = atom[i][j].z_i; q < z1;q++) // calculate ionization rate for z_i<i<z1. DON'T FORGET give current n
				{
					atom[i][j].W_cal(A[q], B[q], fil, n(q+1, k(data_n.ion_pot[q]))); // for direct current; give function n z, which equal q+1, but this equation does't >z1!!!!!!
					atom[i][j].W_cal_AC(A[q], B[q], fil, n(q + 1, k(data_n.ion_pot[q]))); //for alternating current 
				}
				atom[i][j].P();
				if (atom[i][j].j != 0)
				{
					atom[i][j].charge += atom[i][j].j*data_n.charge_step; //Atom's charge will change if degree ionization doesn't equal zero (j!=0)
					atom[i][j].j = 0;
				};
			}
		};

	delete[] A, B, C; //clearing memory from dinamic massiv
}


double C_res(int l, double n)
{
	if (n < 1) { n = 1; }; // not sure 
	return pow(2 * M_PI*n, (-1))*pow((2 * M_E / n), 2 * n); // for l*<<n* e - what is it? electron's mass or exp?
	//return pow(2, 2 * n) / (n*tgamma(n + l + 1)*tgamma(n - l)); //universal
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

//function for initialization of struct (next ficha)
/*
in_data data_he;
void setupStr_l(int value[], int n)
{
for (int i=0; i<n;i++)
data_he.ion_l[i]=value[i];
}
void setupStr_num(int value[], int n)
{
for (int i=0; i<n;i++)
data_he.ion_num_of_el[i]=value[i];
}
void setupStr_pot(double value[], int n)
{
for (int i=0; i<n;i++)
data_he.ion_pot[i]=value[i];*/
