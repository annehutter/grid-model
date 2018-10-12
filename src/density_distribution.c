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
#include "self_shielding.h"
#include "recombination.h"
#include "density_distribution.h"

#define SQR(X) ((X) * (X))

double Hubble(confObj_t simParam)
{
	double Hubble = simParam->h*100; //in km s^-1 Mpc^-1
	return Hubble*sqrt(simParam->omega_m*pow(1.+simParam->redshift,3)+simParam->omega_l);
}

double dd_pdf(double x, void * p)
{
	pdf_params_t * params = (pdf_params_t *)p;
	
	double frac = 2./3.;
	double delta_0 = 7.61/(1.+params->redshift);
	return params->amplitude*exp(-SQR(pow(x,-frac)-params->constant)/(2.*SQR(frac*delta_0)))*pow(x,-params->beta);
}

double dd_calc_integral(pdf_params_t params, double upLim)
{
	gsl_function F;
	F.function = &dd_pdf;
	F.params = &params;
	double result, error;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
	
	if(upLim <=0.)
	{
		gsl_integration_qagiu(&F, 0., 1.e-6, 1.e-6, 10000, w, &result, &error);
	}else{
		gsl_integration_qag(&F, 0., upLim, 0., 1.e-6, 10000, 1, w, &result, &error);
	}
	
	gsl_integration_workspace_free(w);
	return result;
}

double dd_pdf_mass(double x, void * p)
{
	pdf_params_t * params = (pdf_params_t *)p;
	
	double frac = 2./3.;
	double delta_0 = 7.61/(1.+params->redshift);
	return x*params->amplitude*exp(-SQR(pow(x,-frac)-params->constant)/(2.*SQR(frac*delta_0)))*pow(x,-params->beta);
}

double dd_calc_integral_mass(pdf_params_t params, double upLim)
{
	gsl_function F;
	F.function = &dd_pdf_mass;
	F.params = &params;
	double result, error;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
	
	if(upLim <=0.)
	{
		gsl_integration_qagiu(&F, 0., 1.e-6, 1.e-6, 10000, w, &result, &error);
	}else{
		gsl_integration_qag(&F, 0., upLim, 0., 1.e-6, 10000, 1, w, &result, &error);
	}
	
	gsl_integration_workspace_free(w);
	return result;
}

void dd_set_norm_pdf(pdf_params_t * params, double redshift)
{
	double result, result_mass;
	double old_result, old_result_mass;
	double old_amplitude, old_constant, old_beta;
	double old_chi, chi;
	
	int rand_amp, rand_const, rand_beta;
	double precision = 1.e-3;
	
	int counter;
	
	params->redshift = redshift;
	params->amplitude = 0.5;
	params->constant = 0.2;
	params->beta = 2.5;
	
	// compute full integrate
	result = dd_calc_integral(*params, 0.);
	result_mass = dd_calc_integral_mass(*params, 0.);
	chi = (result-1.)*(result-1.)+(result_mass-1.)*(result_mass-1.);

	// adapt struct, so integral is 1
	counter = 0;
	while((fabs(result_mass-1.)>=precision || fabs(result-1.)>=precision) && counter<30000)
	{
		rand_amp = 2*(rand()%2)-1;
		rand_const = 2*(rand()%2)-1;
        rand_beta = 2*(rand()%2)-1;
		
		old_amplitude = params->amplitude;
		old_constant = params->constant;
        old_beta = params->beta;
		old_result = result;
		old_result_mass = result_mass;
		old_chi = chi;
		
		params->amplitude += rand_amp*precision*params->amplitude;
		params->constant += rand_const*precision*params->constant;
        if((1.+rand_beta*precision)*params->beta<=2.5 && (1.+rand_beta*precision)*params->beta>1.)
        {
            params->beta += rand_beta*precision*params->beta;
        }
		result_mass = dd_calc_integral_mass(*params, 0.);
		result = dd_calc_integral(*params, 0.);
		chi = (result-1.)*(result-1.)+(result_mass-1.)*(result_mass-1.);
		
		if(old_chi <= chi)
		{
			params->amplitude = old_amplitude;
			params->constant = old_constant;
            params->beta = old_beta;
			result = old_result;
			result_mass = old_result_mass;
			chi = old_chi;
		}
		counter++;
	}
	printf("counter = %d\n", counter);
	
	printf("%e:\t %e\t%e\t A = %e\t C = %e\t beta = %e\n", redshift, result, result_mass, params->amplitude, params->constant, params->beta);
}

double dd_frac_densSS(double densSS, confObj_t simParam)
{
	double result, result2;
// 	double result_mass;
// 	double old_result, old_result_mass;
// 	double old_amplitude, old_constant;
// 	double old_chi, chi;
	
// 	int rand_amp, rand_const;
// 	double precision = 1.e-4;
	
	// construct struct
	pdf_params_t *params;
	params = malloc(sizeof(pdf_params_t));
	if(params == NULL)
	{
		fprintf(stderr, "params in frac_densSS (density_distribution.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	params->redshift = simParam->redshift;
	params->amplitude = amplitude_norm_pdf(simParam->redshift);
	params->constant = constant_norm_pdf(simParam->redshift);
	params->beta = 2.5;
	
// 	// compute full integrate
// 	result = calc_integral(*params, 0.);
// 	result_mass = calc_integral_mass(*params, 0.);
// 	chi = (result-1.)*(result-1.)+(result_mass-1.)*(result_mass-1.);

// 	// adapt struct, so integral is 1
// 	while(fabs(result_mass-1.)>=precision || fabs(result-1.)>=precision)
// 	{
// 		rand_amp = 2*(rand()%2)-1;
// 		rand_const = 2*(rand()%2)-1;
// 		
// 		old_amplitude = params->amplitude;
// 		old_constant = params->constant;
// 		old_result = result;
// 		old_result_mass = result_mass;
// 		old_chi = chi;
// 		
// 		params->amplitude += rand_amp*precision*params->amplitude;
// 		params->constant += rand_const*precision*params->constant;
// 		result_mass = calc_integral_mass(*params, 0.);
// 		result = calc_integral(*params, 0.);
// 		chi = (result-1.)*(result-1.)+(result_mass-1.)*(result_mass-1.);
// 		
// 		if(old_chi <= chi)
// 		{
// 			params->amplitude = old_amplitude;
// 			params->constant = old_constant;
// 			result = old_result;
// 			result_mass = old_result_mass;
// 			chi = old_chi;
// 		}
// 	}
	


	// compute partial interval & derive fraction
	result = dd_calc_integral(*params, densSS);
	
	result2 = dd_calc_integral(*params, 0.);
	
	free(params);
	
	return result/result2;
}

double dd_calc_mfp(confObj_t simParam, double photHI_bg, double temperature, double redshift)
{
	double lambda_0 = 60./Hubble(simParam);	//in Mpc
	double densSS;
	double mfp;
	
	densSS = ss_calc_densSS(simParam, photHI_bg, temperature, redshift);
	
	mfp = lambda_0*pow(1.-dd_frac_densSS(densSS, simParam),-2./3.); 
	
	return mfp;
}

void set_mfp_Miralda2000(confObj_t simParam)
{
	simParam->mfp = dd_calc_mfp(simParam, simParam->photHI_bg, 1.e4, simParam->redshift);
}

