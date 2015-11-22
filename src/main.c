#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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
#include "sources_to_grid.h"
#include "fraction_q.h"
#include "filtering.h"

int main (int argc, /*const*/ char * argv[]) { 
	int size = 1;
	int myRank = 0;

	char iniFile[1000];
	confObj_t simParam;
	grid_t *grid;
	
	int num_sources;
	source_t *sources;
	
	double t1, t2;
	
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
	
	
	//read paramter file
	simParam = readConfObj(iniFile);
		
	//read files (allocate grid)
	grid = initGrid();
	if(myRank==0) printf("reading files to grid... ",myRank);
	read_files_to_grid(grid, simParam);
	if(myRank==0) printf("done\n");
	
	//read source files
	if(myRank==0) printf("reading sources file... ");
	num_sources = count_sources(simParam->sources_file);
	sources = allocate_source_list(num_sources);
	read_sources(simParam->sources_file, num_sources, sources);
	if(myRank==0) printf("done\n");

	//map sources to grid
	if(myRank==0) printf("mapping sources to grid... ");
	map_nion_to_grid(myRank, grid, num_sources, sources);
	if(myRank==0) printf("done\n");
	
	//compute fraction Q
	if(myRank==0) printf("computing relation between number of ionizing photons and absorptions... ");
	compute_Q(grid, simParam);
	if(myRank==0) printf("done\n");
	
	//apply filtering
	if(myRank==0) printf("apply tophat filter routine for ionization field... ");
	compute_ionization_field(grid);
	if(myRank==0) printf("done\n");
	
	//write ionization field to file
	if(myRank==0) printf("writing ionization field to file... ");
	save_to_file_XHII(grid, simParam->out_XHII_file);
	if(myRank==0) printf("done\n");
	
	//deallocate grid
	if(myRank==0) printf("deallocating grid ...");
	deallocate_grid(grid);
	if(myRank==0) printf("done\n");
	
	//deallocate sources
	if(myRank==0) printf("deallocating sources ...");
	deallocate_source_list(sources);
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