#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

double time_from_redshift_flatuniverse(confObj_t simParam, double zmin, double zmax)
{
	double prefactor = 2./(3*H0*sqrt(simParam->omega_l));
	double tmp = sqrt(simParam->omega_l/simParam->omega_m);
	
	return prefactor*(asinh(tmp*pow(1.+zmin, -1.5)) - asinh(tmp*pow(1.+zmax, -1.5)));
}

void compute_Q(grid_t *thisGrid, confObj_t simParam)
{
	int nbins;
	int local_n0;
	double box_size;
	
	double evol_time;
	double z;
	double mean_numdensity;
	double h;
	
	double Nion, Nabs;
	
	nbins = thisGrid->nbins;
	local_n0 = thisGrid->local_n0;
	box_size = thisGrid->box_size;
	
	if(simParam->calc_ion_history == 1)
	{
		evol_time = simParam->evol_time*Myr_s + time_from_redshift_flatuniverse(simParam, simParam->redshift, simParam->redshift_prev_snap);
		printf("\n zstart = %e\t zend = %e\t evol_time = %e\t %e\n", simParam->redshift_prev_snap, simParam->redshift, simParam->evol_time, time_from_redshift_flatuniverse(simParam, simParam->redshift, simParam->redshift_prev_snap)/Myr_s);
		simParam->evol_time = evol_time/Myr_s;
	}else{
		evol_time = simParam->evol_time*Myr_s;
		printf("\nevol_time = %e\n", evol_time/Myr_s);
	}
	z = simParam->redshift;
	if(simParam->default_mean_density == 1){
		mean_numdensity = rho_g_cm/mp_g*simParam->h*simParam->h*simParam->omega_b*(1.+z)*(1.+z)*(1.+z);
	}else{
		mean_numdensity = simParam->mean_density*(1.+z)*(1.+z)*(1.+z);
	}
	h = simParam->h;
	
	const double volume = pow(box_size/(h*(double)nbins*(1.+z))*Mpc_cm,3);
	
	printf("mean_numdensity = %e\n", mean_numdensity);
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				Nion = creal(thisGrid->nion[i*nbins*nbins+j*nbins+k])*evol_time;
				Nabs = creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity*volume*(1.+creal(thisGrid->nrec[i*nbins*nbins+j*nbins+k]));
				thisGrid->frac_Q[i*nbins*nbins+j*nbins+k] = Nion/Nabs+0.*I;
// 				if(Nion>0.) printf("%e\t%e\t%e\n", Nion, Nabs, Nion/Nabs);
			}
		}
	}
}


