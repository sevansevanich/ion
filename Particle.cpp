#include "Particle.h";
#include "math.h";
#include <iostream>
#include <ctime>
#include <cstdlib>

const double el=1.60217657*pow(10,-19);

Particle::Particle (double ch, int z_1) // конструктор используемый для инициализации объекта (частицы) при создании мишине, с занесением только начального заряда
{
	charge=ch;
	z_i=charge/el;
	z1=z_1; //временно нужно!!! лишняя память по сути, но используется в Рунге-Кутты
}

void Particle::W_cal(int z1, int i, double A, double B, double E) //вычисление частоты ионизации для объекта
{
	W[i]=A*pow(E,(1-2*n))*B*exp(-2/(3*B*E));
}

double Particle::f(double w, double n, int i)
{
	double value;
	if (i<=z_i)
		{
			value=-w*n;
		}
	else
		{
			value=W[i-1]*N[i-1]-w*n;
		}
	if (i==z1-1)
		{
			value=w*n;
		}
	return value;
}

void Particle::N_cal(int M) // метод Рунге-Кутты 4 порядка
{
	double k1,k2,k3,k4,h;

	dt=0.1;
	h=dt/(z1-z_i);

	N[z_i]=1; //в момент времени t=0

	for(int i=z_i;i<z1;i++)
	{
		k1=f(W[i],N[i], i);
		k2=f(W[i]+h/2,N[i]+h/2*k1, i);
		k3=f(W[i]+h/2,N[i]+h/2*k2, i);
		k4=f(W[i]+h,N[i]+h*k3, i);

		N[i+1]=N[i]+(h/6)*(k1+2*k2+2*k3+k4);
	}
}

int Particle::P() // вычисление вероятности ионизации частицы
{
	//std::mt19937 gen(time(0));
	//std::uniform_real_distribution<> urd(0, 1);
	//std::cout << urd(gen) << std::endl;

	double p;
	double sum1,sum2;
	while (j<z1)
	{
		srand(time(0));
		p=rand()/double( RAND_MAX );

		for (int i=z_i;i<=j;i++)
			{
				sum2+=N[i];
				if (i==(j-1)){sum1 =sum2;};
			}
		if ((sum1<p) and (p<sum2)) {return j;};
	};
}
