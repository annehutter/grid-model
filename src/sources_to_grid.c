#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_integration.h>
#include <time.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "confObj.h"
#include "grid.h"
#include "sources.h"

/* map number of ionizing photons to grid --------------------------------------------------------*/
void map_nion_to_grid(int myRank, grid_t *thisGrid, int num_sources, source_t *thisSourceList)
{
	int comx, comy, comz;
	int nbins;
	int local_0_start;
	int local_n0;
	
	nbins = thisGrid->nbins;
	local_0_start = thisGrid->local_0_start;
	local_n0 = thisGrid->local_n0;
	
	for(int source=0; source<num_sources; source++)
	{
		comx = (int)(thisSourceList[source].pos[0]*nbins);
		comy = (int)(thisSourceList[source].pos[1]*nbins);
		comz = (int)(thisSourceList[source].pos[2]*nbins);
				
		if(comx==nbins) comx = comx-1;
		if(comy==nbins) comy = comy-1;
		if(comz==nbins) comz = comz-1;
		
		if(comz>=local_0_start && comz<local_0_start+local_n0)
		{
			thisGrid->nion[(comz-local_0_start)*nbins*nbins + comy*nbins + comx] += thisSourceList[source].Nion*thisSourceList[source].fesc+0.*I;
// 			thisGrid->igm_clump[(comz-local_0_start)*nbins*nbins + comy*nbins + comx] = 1.+0.*I;
		}
	}
}