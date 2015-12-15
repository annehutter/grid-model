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
	double box_size;
	
	double evol_time;
	double z;
	double mean_numdensity;
	double h;
	
	nbins = thisGrid->nbins;
	local_n0 = thisGrid->local_n0;
	box_size = thisGrid->box_size;
	
	evol_time = simParam->evol_time*Myr_s;
	z = simParam->redshift;
	mean_numdensity = simParam->mean_density*(1.+z)*(1.+z)*(1.+z);
	h = simParam->h;
	
	const double volume = pow(box_size/(h*(double)nbins*(1.+z))*Mpc_cm,3);
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
// 				thisGrid->nion[i*nbins*nbins+j*nbins+k] = creal(thisGrid->nion[i*nbins*nbins+j*nbins+k])*evol_time/(creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*creal(thisGrid->igm_clump[i*nbins*nbins+j*nbins+k])*mean_numdensity*pow(box_size/(h*(double)nbins*(1+z))*Mpc_cm,3))+0.*I;
							
				thisGrid->nion[i*nbins*nbins+j*nbins+k] = creal(thisGrid->nion[i*nbins*nbins+j*nbins+k])*evol_time/((1.+creal(thisGrid->igm_clump[i*nbins*nbins+j*nbins+k])*creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity*evol_time*recomb_HII)*creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity*volume)+0.*I;
				
// 				printf("%e\t%e\t%e\t%e\n",creal(thisGrid->igm_clump[i*nbins*nbins+j*nbins+k]),creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity,evol_time,creal(thisGrid->igm_clump[i*nbins*nbins+j*nbins+k])*creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity*evol_time*recomb_HII);
			}
		}
	}
}


