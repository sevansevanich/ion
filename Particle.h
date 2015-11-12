#ifndef PARTICLE_H
#define PARTICLE_H
#include "math.h"

const double el=1.60217657*pow(10,-19);

class Particle

{
	protected:
		void N_cal(int M);
		double f(double,double,int);
		double W[]={}, N[]={};
		double dt;
	public:
		Particle();
		Particle (double, int);
		void W_cal(int, int, double, double, double);
		int P();
		int z_i, j, z1;
		double n, charge;
};

#endif
