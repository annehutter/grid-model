#ifndef MEAN_FREE_PATH_H
#define MEAN_FREE_PATH_H
#endif

typedef struct
{
	double prefactor;
	double prefactor_z;
	double factor;
	double densSS;
	
	double nH;
	double f;
	double dcell;
	double redshift;
} mfp_integrand_t;

double mfp_integrand(double x, void * p);
double calc_mfp_integral(mfp_integrand_t *params);
double calc_local_mfp(confObj_t simParam, double dens, double photHI, double temperature, double redshift);
void compute_web_mfp(grid_t *thisGrid, confObj_t simParam);