#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
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
#include "sources_to_grid.h"
#include "fraction_q.h"
#include "filtering.h"
#include "self_shielding.h"

#include "density_distribution.h"
#include "recombination.h"

void print_mean_photHI(grid_t *thisGrid, confObj_t simParam)
{
  	int nbins;
	int local_n0;
	
	double sum_clump = 0.;
	double sum_XHII = 0.;
	
	double mean_clump;
	double mean_XHII;
	double mean_photHI;
	
	double redshift = simParam->redshift;
	double mean_density = simParam->mean_density*(1.+redshift)*(1.+redshift)*(1.+redshift);
	
	nbins = thisGrid->nbins;
	local_n0 = thisGrid->local_n0;
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				sum_clump += thisGrid->igm_clump[i*nbins*nbins+j*nbins+k];
				sum_XHII += thisGrid->XHII[i*nbins*nbins+j*nbins+k];
			}
		}
	}
	
	mean_clump = sum_clump;
	mean_XHII = sum_XHII;
	
#ifdef __MPI
	MPI_Allreduce(&sum_clump, &mean_clump, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&sum_XHII, &mean_XHII, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	mean_clump = mean_clump/(nbins*nbins*nbins);
	mean_XHII = mean_XHII/(nbins*nbins*nbins);
	
	mean_photHI = mean_XHII*mean_XHII*mean_density*mean_clump*recomb_HII/(1.-mean_XHII);
	printf("dens = %e\t clump = %e\t XHII = %e\t photHI = %e\n",mean_density, mean_clump, mean_XHII, mean_photHI);
}

int main (int argc, /*const*/ char * argv[]) { 
	int size = 1;
	int myRank = 0;

	char iniFile[1000];
	confObj_t simParam;
	grid_t *grid;
	
	sourcelist_t *sourcelist;
	
	integral_table_t *integralTable;
	
	double t1, t2;
	
	double mean_photHI;
	
#ifdef __MPI
	MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &size); 
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank); 
	
	t1 = MPI_Wtime();
#else
	t1 = time(NULL);
#endif
	
	//parse command line arguments and be nice to user
	if (argc != 2) {
		printf("cifog: (C)  - Use at own risk...\n");
		printf("USAGE:\n");
		printf("cifog iniFile\n");
		
		exit(EXIT_FAILURE);
	} else {
		strcpy(iniFile, argv[1]);
	}
	
	//-------------------------------------------------------------------------------
	// reading input files and prepare grid
	//-------------------------------------------------------------------------------
	
	//read paramter file
	simParam = readConfObj(iniFile);
	
	if(simParam->generate_recomb_tables == 1)
	{
		if(myRank==0) printf("generating recombination tables...\n");
		compute_table_norm_pdf(3., 20., 0.01, myRank, size, "norm_pdf.dat");
		compute_table_dens(0.8, 1.2, 0.01, -10., -5., 0.01, myRank, size, "density_integral.dat");
		compute_table_redshift(-2., 2., 0.001, 6., 20., 0.01, simParam, myRank, size, "redshift_integral.dat");
		if(myRank==0) printf("done.\n");
		
#ifdef __MPI
		t2 = MPI_Wtime();
		printf("Execution took %f s\n", t2-t1);
		MPI_Finalize();
#else
		t2 = time(NULL);
		printf("Execution took %f s\n", t2-t1);
#endif
		exit(0);
	}
		
	//read files (allocate grid)
	grid = initGrid();
	if(myRank==0) printf("reading files to grid... ");
	read_files_to_grid(grid, simParam);
	if(myRank==0) printf("done\n");
	
	//read source files (allocate sources)
	if(myRank==0) printf("reading sources file... ");
	sourcelist = read_sources(simParam->sources_file);
	if(myRank==0) printf("done\n");

	//map sources to grid
	if(myRank==0) printf("mapping sources to grid... ");
	map_nion_to_grid(grid, sourcelist);
	if(myRank==0) printf("done\n");
	
	//------------------------------------------------------------------------------
	// compute web model
	//------------------------------------------------------------------------------
	
	if(simParam->use_web_model == 1)
	{
		//set mean free paths
		//this mean free path is an overestimate at high redshifts, becomes correct at z~6
		if(myRank==0) printf("\n\ncompute mean free path... ");
		printf("mean free path at z=%e is %e Mpc... ", simParam->redshift, calc_mfp(simParam));
		if(myRank==0) printf("done\n");
	  
		//set photoionization rate to background value
		if(myRank==0) printf("\n\nsetting photoionization rate to background value... ");
		set_value_to_photoionization_field(grid, simParam);
		if(myRank==0) printf("done\n");
		
// 		if(myRank==0) printf("compute mean photoionization rate... ");
// 		mean_photHI = get_mean_photHI(grid, simParam);
// 		if(myRank==0) printf("done\n");
		
// 		if(simParam->compute_photHIfield == 1)
// 		{
			//write photoionization rate field to file
			if(myRank==0) printf("writing photoionization field to file... ");
			save_to_file_photHI(grid, simParam->out_photHI_file);
			if(myRank==0) printf("done\n");
// 		}
		
		//apply web model
		if(myRank==0) printf("\n\napply web model... ");
		compute_web_ionfraction(grid, simParam);
		if(myRank==0) printf("done\n");
		
		//compute number of recombinations
		if(myRank==0) printf("\n\ncompute number of recombinations... ");
		integralTable = initIntegralTable(6., 20., 0.01, -4., 2., 0.1, -4., 4, 0.1);
		compute_number_recombinations(grid, simParam, "nrec_values_batch_z6_20_0.01_f-4_2_0.1_d-4_4_0.1.dat", integralTable);
		free(integralTable);
		if(myRank==0) printf("done\n");
		
		//compute mean free paths
	}

	//--------------------------------------------------------------------------------
	// apply tophat filter
	//--------------------------------------------------------------------------------
	
	//compute fraction Q
	if(myRank==0) printf("\n\ncomputing relation between number of ionizing photons and absorptions... ");
	compute_Q(grid, simParam);
	if(myRank==0) printf("done\n");
	
	//apply filtering
	if(myRank==0) printf("apply tophat filter routine for ionization field... ");
	compute_ionization_field(grid);
	if(myRank==0) printf("done\n");
	
	//--------------------------------------------------------------------------------
	// writing data to files
	//--------------------------------------------------------------------------------
	
	//write ionization field to file
	if(myRank==0) printf("writing ionization field to file... ");
	save_to_file_XHII(grid, simParam->out_XHII_file);
	if(myRank==0) printf("done\n");
	
	//--------------------------------------------------------------------------------
	// deallocating grids
	//--------------------------------------------------------------------------------
	
	//deallocate grid
	if(myRank==0) printf("deallocating grid ...");
	deallocate_grid(grid);
	if(myRank==0) printf("done\n");
	
	//deallocate sources
	if(myRank==0) printf("deallocating sources ...");
	deallocate_sourcelist(sourcelist);
	if(myRank==0) printf("done\n");
	
	if(myRank==0) printf("Finished\n");
#ifdef __MPI
	t2 = MPI_Wtime();
	printf("Execution took %f s\n", t2-t1);
	MPI_Finalize();
#else
	t2 = time(NULL);
	printf("Execution took %f s\n", t2-t1);
#endif
	
	return 0;
}