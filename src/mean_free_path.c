#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <gsl/gsl_integration.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "phys_const.h"
#include "confObj.h"
#include "grid.h"
#include "sources.h"

#include "density_distribution.h"

#define SQR(X) ((X) * (X))

double calc_local_mfp(pdf_params_t pdf_params, double dens, double XHI, double temperature, double redshift)
{
	const double omega_b = simParam->omega_b;
	const double omega_m = simParam->omega_m;
	const double h = simParam->h;
	const double Y = simParam->Y;
	const double fg = omega_b/omega_m;
	const double mu = 4./(8.-5.*Y);
	double rho;
	
	if(simParam->default_mean_density == 1){
		rho = rho_g_cm;
	}else{
		rho = simParam->mean_density;
	}
	
	tmp = omega_b*H0*clight_cm*SQR(mu)*CUB(mp_g)*G;
	tmp2 = 4.*CUB(M_PI)*SQR(gamma_gas)*SQR(boltzman_cgs)*fg
	
	return tmp/(tmp2*rho*CUB(1.+redshift)*temperature)*XHI*dens*calc_integral(pdf_params, dens);
}