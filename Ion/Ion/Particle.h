#pragma once
#ifndef PARTICLE_H
#define PARTICLE_H
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <iostream>

const double el = 1.60217657*pow(10, -19);
const int kol = 50; // define target's size
const double duration_imp = 40.0; // impulse duration
const double I_h = 13.598434005136; //ionization energy eV
const double E_a = 0.514224*pow(10, 12); // atomic field TV/m
const double w_a = 4.134*pow(10, 16); //atomic frequency
const double m_e = 9.10938291*pow(10, -31); //electon's mass
const double c = 299792458.0; // m/s
const double epsilon0 = 8.8542*pow(10,-12);

class Particle

{
protected:
	std::vector<double> W; // vector uses for dinamic increase memory size 
	std::vector<double> N;
	std::vector<double> W_AC;

	void N_cal_ru();

	void N_cal_dif();
	double P(int, int);
	double P_max(int, int);
public:
	int z_i, j, z1;
	double charge, dt;
	Particle();
	Particle(double, int, int);

	void W_cal(double, double, double, double, double);
	void W_cal_AC(double, double, double, double, double);
	void P();
	void clear();
};

#endif