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
#include "self_shielding.h"
#include "recombination.h"
#include "mean_free_path.h"

#define SQR(X) ((X) * (X))
#define CUB(X) ((X) * (X) * (X))

double pdf_mfp(double x, void * p)
{
	pdf_params_t * params = (pdf_params_t *)p;
	
	double frac = 2./3.;
	double delta_0 = 7.61/(1.+params->redshift);
	return params->amplitude*exp(-SQR(pow(x,-frac)-params->constant)/(2.*SQR(frac*delta_0)))*pow(x,0.5-params->beta);
}

double calc_integral_mfp(pdf_params_t params, double upLim)
{
	gsl_function F;
	F.function = &pdf_mfp;
	F.params = &params;
	double result, error;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
	
	if(upLim <=0.)
	{
		gsl_integration_qagiu(&F, 0., 1.e-9, 1.e-9, 10000, w, &result, &error);
	}else{
		gsl_integration_qag(&F, 0., upLim, 0., 1.e-9, 10000, 1, w, &result, &error);
	}
	
	gsl_integration_workspace_free(w);
	return result;
}

double calc_local_mfp(confObj_t simParam, double dens, double photHI, double temperature, double redshift)
{
	const double omega_b = simParam->omega_b;
	const double omega_m = simParam->omega_m;
	const double Y = simParam->Y;
	double fg = omega_b/omega_m;
	double densSS, modPhotHI;
	double factor, mfp;

	pdf_params_t *params;
	
	params = malloc(sizeof(pdf_params_t));
	if(params == NULL)
	{
		fprintf(stderr, "params in calc_local_mfp (mean_free_path.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	redshift = (1.+redshift)*pow(dens, 1./3.)-1.;
	params->redshift = redshift;
	params->amplitude = amplitude_norm_pdf(redshift);
	params->constant = constant_norm_pdf(redshift);
	params->beta = 2.5;
	
	densSS = ss_calc_densSS(simParam, photHI, temperature, redshift);
// 	printf("densSS = %e\t photHI = %e\t temp = %e\t z = %e\n", densSS, photHI, temperature, redshift);
	modPhotHI = ss_calc_modPhotHI(dens, densSS);
	
	factor = pow(3./(8.*M_PI)*omega_b,0.5)*pow((M_PI*gamma_gas*boltzman_cgs*(8.-5.*Y)*fg*temperature)/(4.*mp_g),-0.5)*pow(1.+redshift,-1.5)*(1.-modPhotHI);
// 	printf("factor = %e\n", factor);
	mfp = 1./(simParam->h*1.e7*CUB(1.+redshift)*factor*calc_integral_mfp(*params,0.));
	
	free(params);
	
	return mfp;
}

void compute_mfp(grid_t *thisGrid, confObj_t simParam)
{
	simParam->mfp = calc_local_mfp(simParam, 1., thisGrid->mean_photHI, 1.e4, simParam->redshift);
	printf("mean free path = %e\n", simParam->mfp);
}

void compute_web_mfp(grid_t *thisGrid, confObj_t simParam)
{
	int nbins;
	ptrdiff_t local_n0;
	int cell;
	double sum=0.;
	double mfp;
	
	double photHI;
	double dens;
	double temperature;
	double redshift;
	
	nbins = thisGrid->nbins;
	local_n0 = thisGrid->local_n0;
	
	temperature = 1.e4;
	redshift = simParam->redshift;

	for(int comz=0; comz<local_n0; comz++)
	{
		for(int comy=0; comy<nbins; comy++)
		{
			for(int comx=0; comx<nbins; comx++)
			{
				cell = comz*nbins*nbins + comy*nbins + comx;
				//compute photHI fluctuations (\delta_{photIon})
				photHI = creal(thisGrid->photHI[cell]);
				dens = creal(thisGrid->igm_density[cell]);

				mfp = calc_local_mfp(simParam, dens, photHI, temperature, redshift);
				sum = sum + mfp;
// 				printf("mfp = %e\n", mfp);
			}
		}
	}
	sum = sum/(local_n0*nbins*nbins);
	printf("mfp = %e\n", sum);
}
