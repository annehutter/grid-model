#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>	//included because M_PI is not defined in <math.h>
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

#include "redshift_tools.h"

double time_from_redshift_flatuniverse(confObj_t simParam, double zmin, double zmax)
{
	double prefactor = 2.*Mpc_cm/(3*simParam->h*1.e7*sqrt(simParam->omega_l));
	double tmp = sqrt(simParam->omega_l/simParam->omega_m);
	
	return prefactor*(asinh(tmp*pow(1.+zmin, -1.5)) - asinh(tmp*pow(1.+zmax, -1.5)));
}
