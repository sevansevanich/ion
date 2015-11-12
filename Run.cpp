#include "Particle.h"

int factorial (int);
double C_res(int, double);
double B_res(double);
double n(int, double);
double k(double);

const int kol=100;
const double I_h=13.598;
const double E_a=1.7*pow(10,7);
const double w_a=4.13*pow(10,16);

struct //  структура входных данных о веществе
{
	double charge=0;
	double charge_step=el;
	int number_of_charge_states = 3;
	int ion_l[2]={0,0};
	int ion_num_of_el[2]={2,1};
	double ion_pot[2]={24.47,54.4};
} data_he;

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
		data_he.ion_pot[i]=value[i];
}*/

int main()
{

	int z0,z1;
	double E;

	z0=data_he.charge/data_he.charge_step;
	z1=z0+data_he.number_of_charge_states-1;

	double A [z1], B [z1], C [z1];

	for (int i=z0;i<=z1;i++) // табулирование коэффициентов для расчета частоты ионизации
		{
			C[i]=C_res(data_he.ion_l[i],n(i,k(data_he.ion_pot[i])));
			B[i]=B_res(data_he.ion_pot[i]);
			A[i]=w_a*(k(pow(data_he.ion_pot[i],2))*factorial(2*data_he.ion_l[i]-1)*pow(B[i],(1-2*n(i,k(data_he.ion_pot[i]))))*pow(C[i],2));
		};

	Particle atom[kol][kol]; // создание мишини

	for (int i=0;i<=kol;i++) // инициализация мишени с начальным зарядом из данных структуры
		for (int j=0;j<=kol;j++)
		{
			atom[i][j].charge=data_he.charge;
			atom[i][j].z1=z1;
		};

	for (int i=0;i<=kol;i++)
		for (int j=0;j<=kol;j++)
		{
			if (atom[i][j].z_i<z1)
			{
				for (int k=atom[i][j].z_i; k<z1;k++) // вычисление частоты ионизации атомы для z_i<i<z1
					atom[i][j].W_cal(z1,k,A[k],B[k],E);
				int res=atom[i][j].P();
				if (res>0) //если нашлось такое j>0, значит частица ионизовалась и изменяем её заряд
					atom[i][j].charge+=atom[i][j].P()*data_he.charge_step; // плюсуем ли????????????????
			}
		};

}

int factorial (int n)
{
  return (n < 2) ? 1 : n * factorial (n - 1);
}
double C_res(int l, double n)
{
	return pow(2,(2*n-2))/(n*gamma(n+l+1)*gamma(n-l));
}
double B_res(double I)
{
	return pow((pow(k(I),3)*E_a),-1);
}
double n(int z, double k)
{
	return z/k;
}
double k(double I)
{
	return sqrt(I/I_h);
}
