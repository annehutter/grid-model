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

#include "utils.h"

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
#include "mean_free_path.h"

#include "input_redshifts.h"
#include "input_grid.h"

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
// 				sum_clump += thisGrid->igm_clump[i*nbins*nbins+j*nbins+k];
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
#ifdef __MPI
	int size = 1;
#endif
	int myRank = 0;

	char iniFile[1000];
	confObj_t simParam;
	
	double *redshift_list;
	
	grid_t *grid = NULL;
	
	sourcelist_t *sourcelist = NULL;
	
	integral_table_t *integralTable = NULL;
	
	double t1, t2;
	
	double zstart = 0., zend = 0., delta_redshift = 0.;
	int snap = -1, num_cycles;
	
	char XHIIFile[1000], cycle_string[8];
	
#ifdef __MPI
	MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &size); 
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank); 
	
	t1 = MPI_Wtime();
	
	fftw_mpi_init();
#else
	t1 = time(NULL);
	
	fftw_init();
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
	
	if(simParam->calc_ion_history == 1){
		num_cycles = simParam->num_snapshots;
	}else{
		num_cycles = 1;
	}
	
	printf("densSS = %e\n", ss_calc_densSS(simParam, 1.e-13, 1.e4, 6.));
	printf("densSS = %e\n", ss_calc_densSS(simParam, 5.1e-11, 1.e4, 14.75));
	printf("densSS = %e\n", ss_calc_densSS(simParam, 1.e-12, 1.e4, 9.));
	printf("densSS = %e\n", ss_calc_densSS(simParam, 1.e-12, 1.e4, 7.)*simParam->omega_b*simParam->h*simParam->h*rho_g_cm/mp_g*8.*8.*8./(1.-simParam->Y));
	printf("z = 6: mfp = %e\n", calc_local_mfp(simParam, 1., 0.5e-12, 1.e4, 6.));
	printf("z = 6: mfp(M2000) = %e\n", dd_calc_mfp(simParam, 0.5e-12, 1.e4, 6.));
	printf("z = 7: mfp = %e\n", calc_local_mfp(simParam, 1., 0.5e-12, 1.e4, 7.));
	printf("z = 7: mfp(M2000) = %e\n", dd_calc_mfp(simParam, 0.5e-12, 1.e4, 7.));
	printf("z = 8: mfp = %e\n", calc_local_mfp(simParam, 1., 0.5e-12, 1.e4, 8.));
	printf("z = 8: mfp(M2000) = %e\n", dd_calc_mfp(simParam, 0.5e-12, 1.e4, 8.));
	printf("z = 9: mfp = %e\n", calc_local_mfp(simParam, 1., 0.5e-12, 1.e4, 9.));
	printf("z = 9: mfp(M2000) = %e\n", dd_calc_mfp(simParam, 0.5e-12, 1.e4, 9.));
	printf("z = 10: mfp = %e\n", calc_local_mfp(simParam, 1., 0.5e-12, 1.e4, 10.));
	printf("z = 10: mfp(M2000) = %e\n", dd_calc_mfp(simParam, 0.5e-12, 1.e4, 10.));
	printf("z = 14.75: mfp = %e\n", calc_local_mfp(simParam, 1., 0.5e-12, 1.e4, 14.75));
	printf("z = 14.75: mfp(M2000) = %e\n", dd_calc_mfp(simParam, 0.5e-12, 1.e4, 14.75));

	//read redshift files with outputs
	redshift_list = NULL;
	if(myRank==0) printf("reading redshift list of files and outputs... ");
	redshift_list = read_redshift_list(simParam->redshift_file, num_cycles);
// 	for(int i=0; i<simParam->num_snapshots; i++) printf("z = %e\t files exist %e\n", redshift_list[2*i], redshift_list[2*i + 1]);
	if(myRank==0) printf("done\n");
	
	//read files (allocate grid)
	grid = initGrid();
	if(myRank==0) printf("reading files to grid... ");
	read_files_to_grid(grid, simParam);
	if(myRank==0) printf("done\n");
	
	if(simParam->calc_ion_history == 1)
	{
		zstart = simParam->redshift_prev_snap;
		zend = simParam->redshift;
		
		if(redshift_list == NULL)
		{
			simParam->redshift_prev_snap = zstart;
			delta_redshift = (zstart-zend)/(double)num_cycles;
			simParam->redshift = zstart - delta_redshift;
		}
	}
	
	for(int cycle=0; cycle<num_cycles - 1; cycle++)
	{
		if(redshift_list != NULL)
		{
			simParam->redshift_prev_snap = redshift_list[2*cycle];
			simParam->redshift = redshift_list[2*(cycle + 1)];
			delta_redshift = redshift_list[2*cycle] - redshift_list[2*(cycle + 1)];
			
			if((redshift_list[2*cycle+1] == 1) || (cycle == 0)) snap++;
		}
		else if(cycle !=0 && redshift_list == NULL)
		{
			simParam->redshift -= delta_redshift;
			simParam->redshift_prev_snap -= delta_redshift;
		}
		
		if(myRank==0) printf("reading sources/nion file for snap = %d... ", snap);
		read_update_nion(simParam, sourcelist, grid, snap);
		if(myRank==0) printf("done\n");
		
		if(myRank==0) printf("reading igm density file... ");
		read_update_igm_density(simParam, grid, snap);
		if(myRank==0) printf("done\n");
	  
		//------------------------------------------------------------------------------
		// compute web model
		//------------------------------------------------------------------------------
		
		if(simParam->use_web_model == 1)
		{
			if(simParam->const_photHI ==1)
			{
				//set photoionization rate to background value
				if(myRank==0) printf("\n\nsetting photoionization rate to background value... ");
				set_value_to_photoionization_field(grid, simParam);
				if(myRank==0) printf("done\n");
			}else{
				if(simParam->calc_mfp == 1)
				{
					//this mean free path is an overestimate at high redshifts, becomes correct at z~6
					if(myRank==0) printf("\n\ncompute mean free path... ");
					set_mfp_Miralda2000(simParam);
					printf("\tmfp = %e\t", simParam->mfp);
					if(myRank==0) printf("done\n");
				}
				
				if(myRank==0) printf("compute mean photoionization rate... ");
				compute_photHI(grid, simParam);
				if(myRank==0) printf("done\n");
			}

			if(simParam->write_photHI_file == 1)
			{
				//write photoionization rate field to file
				if(myRank==0) printf("writing photoionization field to file... ");
				save_to_file_photHI(grid, simParam->out_photHI_file);
				if(myRank==0) printf("done\n");
			}
			
			
			//apply web model
			if(myRank==0) printf("\n\napply web model... ");
			compute_web_ionfraction(grid, simParam);
			if(myRank==0) printf("done\n");
			
			if(simParam->calc_recomb == 1)
			{
				//compute number of recombinations
				if(myRank==0) printf("\n\ncompute number of recombinations... ");
				integralTable = initIntegralTable(simParam->zmin, simParam->zmax, simParam->dz, simParam->fmin, simParam->fmax, simParam->df, simParam->dcellmin, simParam->dcellmax, simParam->ddcell);
				compute_number_recombinations(grid, simParam, simParam->recomb_table, integralTable);
				free(integralTable);
				if(myRank==0) printf("done\n");
			}
			
			if(simParam->calc_mfp == -1)
			{
				//compute mean free paths
				if(myRank==0) printf("\ncompute mean free paths... ");
				compute_web_mfp(grid, simParam);
				if(myRank==0) printf("done\n");
			}
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
		compute_ionization_field(simParam, grid);
		if(myRank==0) printf("done\n");
		
		//write ionization field to file
		for(int i=0; i<100; i++) XHIIFile[i] = '\0';
		strcat(XHIIFile, simParam->out_XHII_file);
		strcat(XHIIFile, "_");
		sprintf(cycle_string,"%02d",cycle); 
		strcat(XHIIFile, cycle_string);
		if(myRank==0) printf("writing ionization field to file %s ... ", XHIIFile);
		save_to_file_XHII(grid, XHIIFile);
		if(myRank==0) printf("done\n");
	}
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
	
	if(myRank==0) printf("deallocating redshift list ...");
	deallocateRedshift_list(redshift_list);
	if(myRank==0) printf("done\n");
	
	confObj_del(&simParam);
	
	if(myRank==0) printf("Finished\n");
#ifdef __MPI
	fftw_mpi_cleanup();
		
	t2 = MPI_Wtime();
	printf("Execution took %f s\n", t2-t1);
	MPI_Finalize();
#else
	fftw_cleanup();
	
	t2 = time(NULL);
	printf("Execution took %f s\n", t2-t1);
#endif
	
	return 0;
}