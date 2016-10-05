#include "Particle.h"
#ifndef DATA_H
#define DATA_H

static struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 2;
	int ion_l[1] = { 0 };
	int ion_n[1] = { 1 };
	int ion_num_of_el[1] = { 1 };
	double ion_pot[1] = { 13.598434005136 };
} data_h;

static struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 3;
	int ion_l[2] = { 0,0 };
	int ion_num_of_el[2] = { 2,1 };
	double ion_pot[2] = { 24.587387936,54.41776311 };
} data_he;

static struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 7;
	int ion_l[6] = { 1,1,0,0,0,0 };
	int ion_n[6] = { 2,2,2,2,1,1 };
	int ion_num_of_el[6] = { 2,1,2,1,2,1 };
	double ion_pot[6] = { 11.26,24.38,47.88,64.49,392.09,489.99 };
} data_c;

static struct
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

static struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 9;
	int ion_l[8] = { 1,1,1,1,0,0,0,0 };
	int ion_n[8] = { 2,2,2,2,2,2,1,1 };
	double ion_pot[8] = { 13.618054, 35.12111, 54.93554, 77.41350, 113.8989, 138.1189, 739.32679, 871.40985 };
	double density = 1.141; //fluid
	double atom_mass = 15.99;
} data_o;

static struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 13;
	int ion_l[12] = { 0,0,1,1,1,1,1,1,0,0,0,0};
	int ion_n[12] = { 3,3,2,2,2,2,2,2,2,2,1,1 };
	double ion_pot[12] = { 7.646235, 15.035267, 80.1436, 109.2654, 141.33, 186.76, 225.02, 265.924, 327.99, 367.489, 1761.80481, 1962.6636 };
	double density = 1.738;
	double atom_mass = 24.307;
} data_mg;

static struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 14;
	int ion_l[13] = { 1,0,0,1,1,1,1,1,1,0,0,0,0 };
	int ion_n[13] = { 3,3,3,2,2,2,2,2,2,2,2,1,1 };
	double ion_pot[13] = { 5.985768, 18.82855, 28.44764, 119.9924, 153.825, 190.49, 241.76, 284.64, 330.21, 398.65, 442.005, 2085.97963, 2304.14 };
	double density = 0.808;
	double atom_mass = 14.01;
} data_al;

static struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 22;
	int ion_l[21] = { 2,2,2,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,0,0 };
	int ion_n[21] = { 3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,1,1 };
	double ion_pot[21] = { 6.56149,12.79977,24.756838,73.4894,91.95,110.68,137.99,158.08,180.03,225.18,249.798,687.36,757.7,833.2,926.5,1008.6,1093.5,1213.1,1287.956,5674.9034,6033.7540 };
	double density = 2.99;
	double atom_mass = 44.95;
} data_sc;

static struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 27;
	int ion_l[26] = { 0,0,2,2,2,2,2,2,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,0,0 };
	int ion_n[26] = { 4, 4, 3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,1,1 };
	double ion_pot[26] = { 7.9024678, 16.19920, 30.651, 54.91, 75.0, 98.985, 124.98, 151.060, 233.6, 262.10, 290.9, 330.8, 361.0, 392.2, 456.2, 489.312, 1262.7, 1357.8, 1460, 1575.6, 1687.0, 1798.4, 1950.4, 2045.759, 8828.1875, 9277.6814 };
	double density = 7.874;
	double atom_mass = 14.01;
} data_fe;

static struct
{
	double charge = 0;
	double charge_step = el;
	int number_of_charge_states = 48;
	int ion_l[47] = { 0,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,0,0,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,0,0 };
	int ion_n[47] = { 5,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,1,1};
	double ion_pot[47] = { 7.576234, 21.4844, 34.83 , 49.0, 65.0, 82.0, 106.0, 125.0, 145.1, 167.0, 271.46, 294.0, 321.0, 347.0, 381.0, 408.43, 469.0, 500.87, 885, 946, 1013, 1082, 1149, 1231, 1308, 1382, 1460, 1535, 1747, 1810.5, 1888.0, 1979, 2077, 2131, 2302, 2371.99, 5558.1, 5753, 5966, 6170, 6551, 6785, 7082, 7271.297, 30097.317, 30965.697 };
	double density = 10.5;
	double atom_mass = 107.8682;
} data_ag;

#endif