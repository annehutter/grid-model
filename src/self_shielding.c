#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>	//included because M_PI is not defined in <math.h>
#include <assert.h>
#include <complex.h>

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

#define SQR(X) ((X) * (X))
#define CUB(X) ((X) * (X) * (X))

double calc_modPhotHI(double densH, double densSS)
{
	double quot;
	
	quot = densH/densSS;
	return 0.98*pow(1.+pow(quot,1.64),-2.28)+0.02*pow(1.+quot,-0.84);
}

double calc_densSS(confObj_t simParam, double photHI, double temperature, double redshift)
{
	double tmp, tmp2, rec_rate;
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
		rho = simParam->mean_density*mp_g;
	}
// 	printf("%e\t%e\t%e\n", photHI, temperature, redshift);
	
	tmp = mu*G/(M_PI*gamma_gas*boltzman_cgs);
	tmp2 = pow(mp_g,5)/fg*pow(1.-Y,-4)/SQR(sigma_HI);
	rec_rate = recomb_HII*pow(temperature*1.e-4,-0.8);	//definition of the HII recombination rate (Fukugita & Kawasaki 1994)
	
	return pow(tmp*tmp2*SQR(photHI)/temperature/SQR(rec_rate),1./3.)/(omega_b*SQR(h)*rho*CUB(1.+redshift));	//checked against Mesinger 2015!
}

double calc_XHII(double dens, double clump, double photHI)
{
	double tmp, tmp2;
	
	tmp = photHI/(dens*clump*recomb_HII);
	tmp2 = 0.5*(tmp + 2. - sqrt((tmp +2.)*(tmp + 2.)-4.));
	if((1.-tmp2)>1.) return 1.;
	else return (1.-tmp2);
}

double calc_photHI_source(source_t *thisSource, double mfp_inv, double boxsize_Mpc, float x, float y, float z)
{
	double r;
	
	const double dx = thisSource->pos[0]-x;
	const double dy = thisSource->pos[1]-y;
	const double dz = thisSource->pos[2]-z;
	
	r = sqrt(dx*dx+dy*dy+dz*dz)*boxsize_Mpc;
	
	return sigma_HI*thisSource->Nion*exp(-r*mfp_inv)/(r*r);
}

void compute_photoionization_field(grid_t *thisGrid, sourcelist_t *thisSourcelist, confObj_t simParam)
{
	int numSources;
	
  	int nbins;
	double nbins_inv, boxsize_cm_inv;
	ptrdiff_t local_0_start, local_n0;
	
	double mfp_inv;
	
	double tmp, sum;
	
	numSources = thisSourcelist->numSources;
	
	nbins = thisGrid->nbins;
	nbins_inv = 1./(double)nbins;
	boxsize_cm_inv = (1.+simParam->redshift)/(simParam->box_size*Mpc_cm);	//inverse boxsize at z in cm

	local_0_start = thisGrid->local_0_start;
	local_n0 = thisGrid->local_n0;
	
	mfp_inv = 1./(simParam->mfp);
	
	sum = 0.;
	for(int comz=0; comz<local_n0; comz++)
	{
		printf("z = %d\n",comz);
	  	const float z = (comz+local_0_start)*nbins_inv;
		for(int comy=0; comy<nbins; comy++)
		{
		  	const float y = comy*nbins_inv;
			for(int comx=0; comx<nbins; comx++)
			{
				const float x = comx*nbins_inv;
				tmp = 0.;
				const source_t *sources = thisSourcelist->source;
				for(int i=0; i<numSources; i++)
				{
					const source_t *thisSource = &sources[i];
					float dx = thisSource->pos[0]-x;
					float dy = thisSource->pos[1]-y;
					float dz = thisSource->pos[2]-z;
					
					if(dx>0.5) dx = dx - 1.0f;
					if(dy>0.5) dy = dy - 1.0f;
					if(dz>0.5) dz = dz - 1.0f;
					if(dx<-0.5) dx = dx + 1.0f;
					if(dy<-0.5) dy = dy + 1.0f;
					if(dz<-0.5) dz = dz + 1.0f;
					
					const float r2 = (dx*dx+dy*dy+dz*dz+0.00001f);
					const float r = sqrtf(r2);
					const float r2_inv = 1.0f/r2;
					
					tmp += sigma_HI*boxsize_cm_inv*thisSource->Nion*thisSource->fesc*boxsize_cm_inv*exp(-r*mfp_inv)*r2_inv;
// 					printf("source: %e\t%e\t%e\n",tmp,r,exp(-r*mfp_inv));
// 				  tmp += calc_photHI_source(thisSourcelist->source[source], mfp_inv, boxsize_Mpc, x, y, z);
				}
				thisGrid->photHI[comz*nbins*nbins + comy*nbins + comx] = tmp + 0.*I;
				sum += tmp;
			}
		}
	}
	thisGrid->mean_photHI = sum;
#ifdef __MPI
	MPI_Allreduce(&sum, &thisGrid->mean_photHI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	thisGrid->mean_photHI = thisGrid->mean_photHI*nbins_inv*nbins_inv*nbins_inv;
}

void construct_photHI_filter(fftw_complex *filter, grid_t *thisGrid, confObj_t simParam)
{
	ptrdiff_t local_n0, local_0_start;
	int nbins;
	
	nbins = thisGrid->nbins;
	local_n0 = thisGrid->local_n0;
	local_0_start = thisGrid->local_0_start;
	
	const double half_nbins = nbins*0.5;
	const double mfp_inv = 1./simParam->mfp;
	const double factor = (thisGrid->box_size/simParam->h/(1.+simParam->redshift))/nbins;
	const double sq_factor = factor*factor;
	
	for(int i=0; i<local_n0; i++)
	{
		const double i_expr = half_nbins-abs(i + local_0_start - half_nbins);		
		const double sq_i_expr = SQR(i_expr);
		for(int j=0; j<nbins; j++)
		{
		  	const double j_expr = half_nbins - abs(j - half_nbins);
			const double sq_j_expr = SQR(j_expr);
			for(int k=0; k<nbins; k++)
			{			
				const double k_expr = half_nbins - abs(k - half_nbins);
				const double sq_k_expr = SQR(k_expr);
			  
				double expr = (sq_i_expr + sq_j_expr + sq_k_expr + 0.25)*sq_factor;
				
				filter[i*nbins*nbins+j*nbins+k] = exp(-sqrt(expr*mfp_inv))/expr + 0.*I;
				
				if(creal(filter[i*nbins*nbins+j*nbins+k])<=0.) printf("%d: %e\n",(int)local_0_start,creal(filter[i*nbins*nbins+j*nbins+k]));
			}
		}
	}
#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void convolve_fft_photHI(grid_t *thisGrid, fftw_complex *filter, fftw_complex *nion_smooth)
{
	int nbins;
	double factor;
	
#ifdef __MPI
	ptrdiff_t alloc_local, local_n0, local_0_start;
	ptrdiff_t n0, n1, n2, local_n1, local_1_start, local_n0x, local_0x_start;
#else
	ptrdiff_t local_n0;
#endif
	
	double mean_photHI;
	
	fftw_complex *nion;
	fftw_complex *nion_ft, *filter_ft;
	fftw_plan plan_nion, plan_filter, plan_back;
	
	nbins = thisGrid->nbins;
	nion = thisGrid->nion;
	local_n0 = thisGrid->local_n0;
	
#ifdef __MPI
	local_0_start = thisGrid->local_0_start;
	
	n0 = nbins;
	n1 = nbins;
	n2 = nbins;
	
	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
	assert(local_0_start == thisGrid->local_0_start);
	assert(local_n0 == thisGrid->local_n0);
	
	nion_ft = fftw_alloc_complex(alloc_local);
	filter_ft = fftw_alloc_complex(alloc_local);
	
	plan_nion = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, nion, nion_ft, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_filter = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, filter, filter_ft, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);
#else 
	nion_ft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	filter_ft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	
	plan_nion = fftw_plan_dft_3d(nbins, nbins, nbins, nion, nion_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_filter = fftw_plan_dft_3d(nbins, nbins, nbins, filter, filter_ft, FFTW_FORWARD, FFTW_ESTIMATE);
#endif
	
	fftw_execute(plan_nion);
	fftw_execute(plan_filter);

#ifdef __MPI
	fftw_mpi_local_size_3d_transposed(n0, n1, n2, MPI_COMM_WORLD, &local_n0x, &local_0x_start, &local_n1, &local_1_start);
#endif
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				nion_ft[i*nbins*nbins+j*nbins+k] = nion_ft[i*nbins*nbins+j*nbins+k]*filter_ft[i*nbins*nbins+j*nbins+k];
			}
		}
	}
	
#ifdef __MPI
	plan_back = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, nion_ft, nion_smooth, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);
#else
	plan_back = fftw_plan_dft_3d(nbins, nbins, nbins, nion_ft, nion_smooth, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
	fftw_execute(plan_back);
	
	factor = 1./(nbins*nbins*nbins);
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				nion_smooth[i*nbins*nbins+j*nbins+k] = factor*nion_smooth[i*nbins*nbins+j*nbins+k];
			}
		}
	}
	
	mean_photHI = creal(nion_smooth[0])*sigma_HI/(Mpc_cm*Mpc_cm);
	
	fftw_destroy_plan(plan_nion);
	fftw_destroy_plan(plan_filter);
	fftw_destroy_plan(plan_back);
	
	fftw_free(nion_ft);
	fftw_free(filter_ft);
	
#ifdef __MPI
	MPI_Bcast(&mean_photHI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
	printf("mean photHI = %e\n", mean_photHI);
	
	thisGrid->mean_photHI = mean_photHI;
}

void compute_photHI(grid_t *thisGrid, confObj_t simParam)
{
#ifdef __MPI
	ptrdiff_t alloc_local, local_n0, local_0_start;
#else
	ptrdiff_t local_n0;
#endif
	int nbins = thisGrid->nbins;
	fftw_complex *filter;
	fftw_complex *nion_smooth;
	
	const double factor=sigma_HI/(SQR(Mpc_cm));

	
#ifdef __MPI
	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
	filter = fftw_alloc_complex(alloc_local);
	nion_smooth = fftw_alloc_complex(alloc_local);
	assert(local_n0 == thisGrid->local_n0);
	assert(local_0_start == thisGrid->local_0_start);
#else
	local_n0 = thisGrid->local_n0;
	filter = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	nion_smooth = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
#endif
	
	//construct exp(-r/mfp)/(r*r) filter
	construct_photHI_filter(filter, thisGrid, simParam);
	
#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				nion_smooth[i*nbins*nbins+j*nbins+k] = creal(thisGrid->nion[i*nbins*nbins+j*nbins+k]) + 0.*I;
			}
		}
	}
	
	//apply filter to Nion field
	convolve_fft_photHI(thisGrid, filter, nion_smooth);
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				thisGrid->photHI[i*nbins*nbins+j*nbins+k] = creal(nion_smooth[i*nbins*nbins+j*nbins+k])*factor + 0.*I;
			}
		}
	}
	
	fftw_free(filter);
	fftw_free(nion_smooth);
}

void set_value_to_photoionization_field(grid_t *thisGrid, confObj_t simParam)
{
	ptrdiff_t local_n0;
	int nbins;
	
	local_n0 = thisGrid->local_n0;
	nbins = thisGrid->nbins;
	
	initialize_grid(thisGrid->photHI, nbins, local_n0, simParam->photHI_bg);
	thisGrid->mean_photHI = simParam->photHI_bg;
}

void compute_web_ionfraction(grid_t *thisGrid, confObj_t simParam)
{
  	int nbins;
	ptrdiff_t local_n0;
	int cell;
	
	double mean_density;
	double photHI;
	double densSS;
	double mod_photHI;
	double redshift;
	double temperature;
	
	nbins = thisGrid->nbins;
	local_n0 = thisGrid->local_n0;
	
	redshift = simParam->redshift;
	temperature = 1.e4;
	mean_density = simParam->mean_density*(1.+redshift)*(1.+redshift)*(1.+redshift);
	
	//compute modified photoionization fraction & compute new ionization fraction in ionized regions

	for(int comz=0; comz<local_n0; comz++)
	{
		for(int comy=0; comy<nbins; comy++)
		{
			for(int comx=0; comx<nbins; comx++)
			{
				cell = comz*nbins*nbins + comy*nbins + comx;
					//compute photHI fluctuations (\delta_{photIon})
					photHI = creal(thisGrid->photHI[cell])*1.e12;
					
					//compute self shielded overdensity
					densSS = calc_densSS(simParam, photHI, temperature, redshift);
					
					//compute modified photHI
					mod_photHI = calc_modPhotHI(creal(thisGrid->igm_density[cell]), densSS);
					
					//compute new XHII
					thisGrid->XHII[cell] = calc_XHII(creal(thisGrid->igm_density[cell])*mean_density, creal(thisGrid->igm_clump[cell]), mod_photHI*creal(thisGrid->photHI[cell])) + 0.*I;
			}
		}
	}
}

