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

#include "self_shielding.h"
#include "mean_free_path.h"

#define SQR(X) ((X) * (X))
#define CUB(X) ((X) * (X) * (X))

double mfp_integrand(double x, void * p)
{
	mfp_integrand_t * params = (mfp_integrand_t *)p;
	double term, term2;
	
	term = exp(-0.019426*SQR(1.+params->redshift)*SQR(pow(x,-2./3.)-1.+6.68589*exp(-0.66*params->redshift))*pow(params->dcell,2./3.))/CUB(x);
	term2 = 1.-2./(1.+pow(1.+4.*params->nH*x*params->f*CUB(1.+params->redshift),0.5));
	
// 	printf("%e\n",params->prefactor*params->prefactor_z*params->factor*term*term2);
	return params->prefactor*params->prefactor_z*params->factor*term*term2*(1.-calc_modPhotHI(x, params->densSS));
}

double calc_mfp_integral(mfp_integrand_t *params)
{
	gsl_function F;
	F.function = &mfp_integrand;
	F.params = (void *)params;
	double result, error;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
	
	gsl_integration_qagiu(&F, 0., 1.e-9, 1.e-9, 10000, w, &result, &error);
	
	gsl_integration_workspace_free(w);

	return result;
}

double calc_local_mfp(confObj_t simParam, double dens, double photHI, double temperature, double redshift)
{
	const double omega_b = simParam->omega_b;
	const double omega_m = simParam->omega_m;
	const double Y = simParam->Y;
	double fg = omega_b/omega_m;
	double nH;
	double f;
	double mfp;
	
	if(simParam->default_mean_density == 1){
		nH = 3.*SQR(H0)/(8.*M_PI*G)/mp_g*omega_b*(1.-Y);
	}else{
		nH = simParam->mean_density*(1.-Y);
	}
	f = recomb_HII/photHI*(1.-0.75*Y)/(1.-Y);
	
	mfp_integrand_t *mfp_integrand;
	
	mfp_integrand = malloc(sizeof(mfp_integrand_t));
	if(mfp_integrand == NULL)
	{
		fprintf(stderr, "mfp_integrand in calc_local_mfp (mean_free_path.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	mfp_integrand->prefactor = 3*clight_cm*H0*pow(CUB(1-Y)*fg*G,0.5)*omega_b/(4.*G*pow(M_PI,1.5)*pow((1.-0.5*Y)*boltzman_cgs*(8.-5.*Y)*gamma_gas,0.5));
	mfp_integrand->prefactor_z = (1.+6.68589*exp(-0.72*redshift))*(0.03+0.053*redshift)*pow(1.+redshift,-4.5);
	mfp_integrand->factor = pow(SQR(f)*CUB(nH)*temperature,-0.5);
	
	mfp_integrand->densSS = calc_densSS(simParam, photHI, temperature, redshift);
	
	mfp_integrand->nH = nH;
	mfp_integrand->f = f;
	mfp_integrand->dcell = dens;
	mfp_integrand->redshift = redshift;

	mfp = clight_cm/(H0*CUB(1.+redshift)*calc_mfp_integral(mfp_integrand)*Mpc_cm);
	
	free(mfp_integrand);
	
	return mfp;
}

void compute_web_mfp(grid_t *thisGrid, confObj_t simParam)
{
	int nbins;
	ptrdiff_t local_n0;
	int cell;
	double sum=0.;
	
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
					
				sum += calc_local_mfp(simParam, dens, photHI, temperature, redshift);
			}
		}
	}
	sum = sum/(local_n0*nbins*nbins);
	printf("mfp = %e\n", sum);
}