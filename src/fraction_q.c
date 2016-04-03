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

void compute_Q(grid_t *thisGrid, confObj_t simParam)
{
	int nbins;
	int local_n0;
// 	int local_0_start;
	double box_size;
	
	double evol_time;
	double z;
	double mean_numdensity;
	double h;
	
	double Nion, Nabs;
	
	nbins = thisGrid->nbins;
	local_n0 = thisGrid->local_n0;
// 	local_0_start = thisGrid->local_0_start;
	box_size = thisGrid->box_size;
	
	evol_time = simParam->evol_time*Myr_s;
	z = simParam->redshift;
	if(simParam->default_mean_density == 1){
		mean_numdensity = rho_g_cm/mp_g*0.2*(1.+z)*(1.+z)*(1.+z)*0.01;
	}else{
		mean_numdensity = simParam->mean_density*(1.+z)*(1.+z)*(1.+z);
	}
	h = simParam->h;
	
	const double volume = pow(box_size/(h*(double)nbins*(1.+z))*Mpc_cm,3);
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
// 				thisGrid->nion[i*nbins*nbins+j*nbins+k] = creal(thisGrid->nion[i*nbins*nbins+j*nbins+k])*evol_time/(creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*creal(thisGrid->igm_clump[i*nbins*nbins+j*nbins+k])*mean_numdensity*pow(box_size/(h*(double)nbins*(1+z))*Mpc_cm,3))+0.*I;
				Nion = creal(thisGrid->nion[i*nbins*nbins+j*nbins+k])*evol_time;
// 				Nabs = ((1.+creal(thisGrid->igm_clump[i*nbins*nbins+j*nbins+k])*creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity*evol_time*recomb_HII)*creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity*volume);
				Nabs = creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity*volume*(1.+creal(thisGrid->nrec[i*nbins*nbins+j*nbins+k]));
// 				if(Nion>0.) printf("%e\t%e\t%e\t%e\t%e\n",Nion,Nabs,volume, mean_numdensity, creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k]));
				thisGrid->nion[i*nbins*nbins+j*nbins+k] = Nion/Nabs+0.*I;
				
// 				printf("%e\t%e\t%e\t%e\n",creal(thisGrid->igm_clump[i*nbins*nbins+j*nbins+k]),creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity,evol_time,creal(thisGrid->igm_clump[i*nbins*nbins+j*nbins+k])*creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity*evol_time*recomb_HII);
			}
		}
	}
	
// 	write_grid_to_file_float(thisGrid->nion, nbins, local_n0, local_0_start, "Q.out");
}


