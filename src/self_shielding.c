#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

double calc_modPhotHI(double densH, double densSS)
{
	double quot;
	
	quot = densH/densSS;
	return 0.98*pow(1.+pow(quot,1.64),-2.28)+0.02*pow(1.+quot,-0.84);
}

double calc_densSS(double photHI, double temperature, double redshift)
{
	return 15.*pow(photHI*10.,2./3.)*pow(temperature*1.e-4,-0.13)*pow((1.+redshift)/7.,-3);
}

double calc_XHII(double dens, double clump, double photHI)
{
	double tmp, tmp2;
	
// 	tmp = 2.*dens*clump*recomb_HII;
// 	tmp2 = photHI/tmp*(sqrt(1.+2.*tmp/photHI)-1.);
	
	tmp = photHI/(dens*clump*recomb_HII);
	tmp2 = 0.5*(tmp + 2. - sqrt((tmp +2.)*(tmp + 2.)-4.));
	printf("XHII = %e\t dens = %e\t clump = %e\t photHI = %e\n", tmp2, dens, clump, photHI);
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
	double nbins_inv, boxsize_cm, boxsize_cm_inv;
	ptrdiff_t local_0_start, local_n0;
	
	double mfp_inv;
	
	double tmp, sum;
	
	numSources = thisSourcelist->numSources;
	
	nbins = thisGrid->nbins;
	nbins_inv = 1./(double)nbins;
	boxsize_cm = simParam->box_size/(1.+simParam->redshift)*Mpc_cm;	//boxsize at z in cm
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
// 	printf("sum = %e\t sum photHI = %e\n", sum, thisGrid->mean_photHI);
#endif
	thisGrid->mean_photHI = thisGrid->mean_photHI*nbins_inv*nbins_inv*nbins_inv;
// 	printf("mean photHI = %e\n", thisGrid->mean_photHI);
}

void set_value_to_photoionization_field(grid_t *thisGrid, confObj_t simParam)
{
	ptrdiff_t local_n0, local_0_start;
	int nbins;
	
	local_n0 = thisGrid->local_n0;
	local_0_start = thisGrid->local_0_start;
	nbins = thisGrid->nbins;
	
	initialize_grid(thisGrid->photHI, nbins, local_n0, local_0_start, simParam->photHI_bg);
	thisGrid->mean_photHI = simParam->photHI_bg;
}

void compute_web_ionfraction(grid_t *thisGrid, confObj_t simParam)
{
  	int nbins;
	ptrdiff_t local_n0;
	int cell;
	
	double mean_density;
	double photHI;
	double mean_photHI_inv;
	double densSS;
	double mod_photHI;
	double redshift;
	double temperature;
	
	nbins = thisGrid->nbins;
	local_n0 = thisGrid->local_n0;
	
	mean_photHI_inv = 1./thisGrid->mean_photHI*1.e-12;
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
				if(creal(thisGrid->XHII[cell])==1.)
				{
					//compute photHI fluctuations (\delta_{photIon})
					photHI = creal(thisGrid->photHI[cell])*1.e12;
					
					//compute self shielded overdensity
					densSS = calc_densSS(photHI, temperature, redshift);
					
					//compute modified photHI
					mod_photHI = calc_modPhotHI(creal(thisGrid->igm_density[cell]), densSS);
// 					printf("%e\t%e\n",densSS, mod_photHI*creal(thisGrid->mean_photHI));
					
					//compute new XHII
					thisGrid->XHII[cell] = calc_XHII(creal(thisGrid->igm_density[cell])*mean_density, creal(thisGrid->igm_clump[cell]), mod_photHI*creal(thisGrid->mean_photHI)) + 0.*I;
				}
			}
		}
	}
}

