#include "Particle.h"
#include <random> // for GSCH
#include <ctime> // fo GSCH
#include <cstdlib>

std::mt19937 gen(time(0));
std::uniform_real_distribution<> urd(0, 1);

Particle::Particle() //consruction by default
{
	z_i = 0;
	j = 0;
	z1 = 0;
	charge = 0;
	dt = 0.2668*pow(10, -14); //laser period for wavelength = 0.8 mkm or this must be the time step in PIC code?
}

Particle::Particle(double ch, int z_1, int n_p)// initialization particle from input struct
{
	charge = ch; //particle's charge
	z_i = charge / el; // quantity ionization
	z1 = z_1; //  temporary arrangement !!!
	j = 0; // amount ionizationing electrons on step, on first it equals zero. reset to zero for next step
	dt = 0.2668*pow(10, -14); //laser period for wavelength = 0.8 mkm or this must be the time step in PIC code? here we increase to -3 degree or else don't correctly work
}

void Particle::W_cal(double A, double B, double E, double n, double keld)//ionization rate ADK for DC
{
	W.insert(W.end(), A * pow(E, -(2 * n - 1))*exp((-2 * B / (3 * E))*(1 - 0.1*keld*keld))); // for vector massive
}

void Particle::W_cal_AC(double A, double B, double E, double n, double keld)//ionization rate ADK for AC
{
	W_AC.insert(W_AC.end(), sqrt((3*E)/(M_PI*B))*A*pow(E, -(2 * n - 1))*exp((-2 * B / (3 * E))*(1 - 0.1*keld*keld)));
}


double Particle::f(double w, double n, int i)//support function
{
	double value;
	if (i <= z_i)
	{
		value = -w*n;
	}
	else
	{
		value = W.at(i - 1) * N.at(i - 1) - w*n;
	}
	if (i == z1 - 1)
	{
		value = w*n;
	}
	return value;
}

void Particle::N_cal_ru() // method of Runge-Kutti 4 degree
{
	double k1, k2, k3, k4, h;

	h = dt / (z1 - z_i);

	N.insert(N.end(), 1);  //initial date t=0

	for (int i = z_i;i<z1;i++)
	{
		k1 = f(W.at(i), N.at(i), i);
		k2 = f(W.at(i) + h / 2, N.at(i) + h / 2 * k1, i);
		k3 = f(W.at(i) + h / 2, N.at(i) + h / 2 * k2, i);
		k4 = f(W.at(i) + h, N.at(i) + h*k3, i);

		N.insert(N.end(), N.at(i) + (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4));
	}
}

void Particle::N_cal_dif() //see when ionizate from z_i=0 to k_max!!!!
{
	double p = 0;
	N.insert(N.end(), 0); // initialization for sum function -1 eq.(11)
	N.at(0) += exp(-W.at(z_i) * dt); // probability not ionize the ion
	for (int k = 1;k <= z1-z_i;k++) // while k (degree ionization at that step) < max charge calculate coef k [1;z1-z_i] (k_max-i)
		{
			N.insert(N.end(), 0); //add to vector new element 
			// certainly read!! -> i=z, z=z_i+1 - particle charge after once ionization?, but we mustn't +1 in function P and P_max because our massive W begin from zero 0 
			if (z_i == z1-1) { p += 1 - exp(-W.at(z1-1)*dt); } // if start state equal k_max-1 (one possible ionization) => k must equal 1
			else {
				if (k == z1-z_i) //check for max value of ionization's ratio
				{
					for (int h = 1;h <= z1 - (z_i + 1);h++) //sum for eq.(8), where P_max is the function in sum. h [1;k_max-i] i=z_i+1 k_max=z1 - max charge
					{
						p += P_max(k, h);
					};
				}
				else
				{
					for (int h = 1;h <= k;h++) //sum for eq.(7), where P is the function in sum. h [1;k_max-i] and k!=k_max!=z1
					{
						p += P(k, h);
					};
				};
			}
			N.at(k)+=p; //N=p in eq.(7)? F from eq.(10) calculate at P
			p = 0;
		}
}

double Particle::P(int k_in, int p)  //probability for i+k degree ionization for N_cal_dif (from eq.(7)) for N_cal_dif
{
	//double k = k_in - 1; // for save vaector's boundary 
	double m(0), res(1), prob_ion_k(0), prob_ion_p(0);
	prob_ion_k = W.at(k_in + z_i);
	prob_ion_p = W.at(p + z_i - 1);
	m = (exp(-prob_ion_k*dt) - exp(-prob_ion_p*dt)) / (1 - prob_ion_k / prob_ion_p);
	for (int j = 1;j <= k_in;j++)
		if (j != p) { res *= 1 - prob_ion_p / W.at(j + z_i - 1); };
	res = 1 / res;
	res *= m;
	return res;
}

double Particle::P_max(int k_in, int p)
{
	double k = k_in - 1;
	double m(0), res(1), prob_ion_k(0), prob_ion_p(0);
	prob_ion_k = W.at(k + z_i);
	prob_ion_p = W.at(p + z_i - 1);
	m = 1+ ((prob_ion_k/prob_ion_p)*exp(-prob_ion_p*dt) - exp(-prob_ion_k*dt)) / (1 - prob_ion_k / prob_ion_p); // difference between P and P_max
	for (int j = 1;j <= k-z_i;j++)
		if (j != p) { res *= 1 - prob_ion_p / W.at(j + z_i - 1); }; // difference between P and P_max
	res = 1 / res;
	res *= m;
	return res;
}

void Particle::P() // void because j - number of new electrons and degree of ionization reside in the object. DON'T FORGET set it zero after create new particle
{
	//N_cal_ru();
	N_cal_dif();

	double sum1 = 0, sum2 = 0;
	int k = -1, exit = 0; // exit - for exit out loop 
	
	while (exit == 0)
	{
		k++;
		if (k > z1-z_i) { exit = 1; } //condition to exit, when k=z1
		else
		{
			double u = urd(gen);
			sum1 = sum2;
			sum2 += N.at(k);
			if ((u > sum1) && (u < sum2)) { j = k; z_i += j; exit = 1; };
		}
	};
}

void Particle::clear()
{
	W.clear();
	W_AC.clear();
	N.clear();
}