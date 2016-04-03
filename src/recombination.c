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
#include "grid.h"
#include "sources.h"

#include "density_distribution.h"
#include "recombination.h"

#define SQR(X) ((X) * (X))

/* two cases for computing recombinations:
 1) only one snapshot available
 2) multiple snapshots available, just compute recombinations for that step and add them on nrec
*/

// double calc_XHI(double dens, double photHI)
// {
// 	double tmp;
// 	
// 	tmp = photHI/(2.*recomb_HII*dens);
// 	return 1. + tmp - sqrt((tmp*tmp)-1.);
// }

void set_nrec_params(nrec_params_t *nrec_params, pdf_params_t *pdf_params, confObj_t simParam, double dens, double photHI, double temp, double redshift)
{
	nrec_params->dens = dens;
	nrec_params->photHI = photHI;
	nrec_params->temp = temp;
	nrec_params->redshift = redshift;
	nrec_params->simParam = simParam;
// 	printf("test: %e\t%e\n", nrec_params->photHI, nrec_params->dens);

// 	set_norm_pdf(pdf_params, simParam);
}

double nrec_rate(double x, void * p)
{
	nrec_params_t * params = (nrec_params_t *)p;
	  
	double frac = 2./3.;
	double delta_0 = 7.61/(1.+params->redshift);
	double nrec;
	
	nrec = 0.5*params->photHI*(sqrt(1.+4.*x*params->dens*recomb_HII/params->photHI)-1.);
	printf("%e\t%e\t%e\t%e\n", params->photHI, params->dens, recomb_HII, x);
	printf("1: %e\t%e\n", params->amplitude*exp(-SQR(pow(x,-frac)-params->constant)/(2.*SQR(frac*delta_0)))*pow(x,-params->beta), nrec);
	
	return params->amplitude*exp(-SQR(pow(x,-frac)-params->constant)/(2.*SQR(frac*delta_0)))*pow(x,-params->beta)*nrec;
}

double calc_nrec_rate(nrec_params_t *params, double redshift)
{
	params->redshift = redshift;
	
	gsl_function F;
	F.function = &nrec_rate;
	F.params = &params;
	double result, error;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
	
	gsl_integration_qagiu(&F, 0., 1.e-9, 1.e-9, 10000, w, &result, &error);
	
	gsl_integration_workspace_free(w);
	printf("2: %e\n", result);
	return result;
}

double calc_nrec_rate_integrand(double x, void * p)
{
	nrec_params_t * params = (nrec_params_t *)p;
	confObj_t simParam = params->simParam;
	
	double tmp = H0*(simParam->omega_b*(1.+x)*(1.+x)*(1.+x)+simParam->omega_l)*(1.+x);
	
	printf("3: %e\n", calc_nrec_rate(params, x)/tmp);

	return calc_nrec_rate(params, x)/tmp;
}

double calc_nrec_integral(nrec_params_t *params, double zstart, double zend)
{
	gsl_function F;
	F.function = &calc_nrec_rate_integrand;
	F.params = &params;
	double result, error;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
	
	gsl_integration_qag(&F, zstart, zend, 0., 1.e-9, 10000, 1, w, &result, &error);
	
	gsl_integration_workspace_free(w);
	
	printf("4: %e\n", result);
	return result;
}

double get_nrec_history(nrec_params_t *nrec_params, double zstart, double redshift)
{
	printf("%e\n", calc_nrec_integral(nrec_params, zstart, redshift));
	return calc_nrec_integral(nrec_params, zstart, redshift);
}

double get_nrec_timestep(nrec_params_t *nrec_params, double zstart, double zend)
{
	return calc_nrec_integral(nrec_params, zstart, zend);
}

void compute_number_recobinations(grid_t *thisGrid, confObj_t simParam)
{
  	int nbins;
	int local_n0;
	double box_size;
	
	double dens;
	double photHI;
	double temp;
	double zstart;
	double redshift;
	
	nrec_params_t *nrec_params;
	pdf_params_t *pdf_params;
	
	nbins = thisGrid->nbins;
	local_n0 = thisGrid->local_n0;
	box_size = thisGrid->box_size;
	
	temp = 1.e4;
	zstart = 10.;
	redshift = simParam->redshift;
	
	pdf_params = malloc(sizeof(pdf_params_t));
	nrec_params = malloc(sizeof(nrec_params_t));
	
	set_norm_pdf(pdf_params, simParam);

	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				dens = creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k]);
				photHI = creal(thisGrid->photHI[i*nbins*nbins+j*nbins+k]);
				set_nrec_params(nrec_params, pdf_params, simParam, dens, photHI, temp, redshift);
				thisGrid->nrec[i*nbins*nbins+j*nbins+k] = 0.+0.*I;// get_nrec_history(nrec_params, zstart, redshift)+0.*I;
			}
		}
	}
	
	free(pdf_params);
	free(nrec_params);
}
