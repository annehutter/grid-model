#include <stdio.h>
#include <math.h>
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
void map_nion_to_grid(grid_t *thisGrid, sourcelist_t *thisSourcelist)
{
	int num_sources;
	int comx, comy, comz;
	int nbins;
	int local_0_start;
	int local_n0;
	
	int i=0;
	
	num_sources = thisSourcelist->numSources;
	
	nbins = thisGrid->nbins;
	local_0_start = thisGrid->local_0_start;
	local_n0 = thisGrid->local_n0;
	
	for(source_t *source=thisSourcelist->source; i<num_sources; i++, source++)
	{
		comx = (int)(source->pos[0]*nbins);
		comy = (int)(source->pos[1]*nbins);
		comz = (int)(source->pos[2]*nbins);
				
		if(comx==nbins) comx = comx-1;
		if(comy==nbins) comy = comy-1;
		if(comz==nbins) comz = comz-1;
		
		if(comz>=local_0_start && comz<local_0_start+local_n0)
		{
			thisGrid->nion[(comz-local_0_start)*nbins*nbins + comy*nbins + comx] += source->Nion*source->fesc+0.*I;
// 			thisGrid->igm_clump[(comz-local_0_start)*nbins*nbins + comy*nbins + comx] = 1.+0.*I;
		}
	}
}