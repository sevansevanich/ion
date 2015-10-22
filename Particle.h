class Particle

{
	protected:
		void W_cal(int, int, double, double, double);
		void N_cal(int M);
		double f(double,double,int);
		double W[]={}, N[]={};
		double dt;
	public:
		Particle (double, int);
		int P();
		int z_i, j, z1;
		double n, charge;
};
