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
#include "photion_background.h"
#include "sources_to_grid.h"
#include "fraction_q.h"
#include "filtering.h"
#include "self_shielding.h"

#include "density_distribution.h"
#include "recombination.h"
#include "mean_free_path.h"

#include "input_redshifts.h"
#include "input_grid.h"
#include "checks.h"

#include "cifog.h"

int main (int argc, /*const*/ char * argv[]) { 
    int size = 1;
    int myRank = 0;

    char iniFile[MAXLENGTH];
    confObj_t simParam;
    
    double *redshift_list = NULL;
    
    grid_t *grid = NULL;
    
    sourcelist_t *sourcelist = NULL;
    
    integral_table_t *integralTable = NULL;
    
    photIonlist_t *photIonBgList = NULL;
    
    double t1, t2;
    
    int num_cycles = 0;
    int restart = 0;
        
#ifdef __MPI
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank); 
    
    t1 = MPI_Wtime();
    
    fftw_mpi_init();
#else
    t1 = time(NULL);
#endif
    
    //parse command line arguments and be nice to user
    if (argc != 2 && argc != 3) {
        printf("cifog: (C)  - Use at own risk...\n");
        printf("USAGE:\n");
        printf("cifog iniFile\n");
        
        printf("argc = %d\n", argc);
        
        exit(EXIT_FAILURE);
    } else {
        if(argc == 2)
        {
            strcpy(iniFile, argv[1]);
        }else{
            if(strcmp(argv[1], "-c") == 0)
            {
                restart = 1;
            }
            strcpy(iniFile, argv[2]);
        }
    }

    cifog_init(iniFile, &simParam, &redshift_list, &grid, &integralTable, &photIonBgList, &num_cycles, restart, myRank);
    
    cifog(simParam, redshift_list, grid, sourcelist, integralTable, photIonBgList, num_cycles, myRank, size);
    
    cifog_deallocate(simParam, redshift_list, grid, integralTable, photIonBgList, myRank);

#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    if(myRank==0) printf("\n********\nFINISHED\n********\n");
    
#ifdef __MPI
    fftw_mpi_cleanup();
        
    t2 = MPI_Wtime();
    if(myRank == 0) printf("\nExecution took %f s\n", t2-t1);
    MPI_Finalize();
#else
    fftw_cleanup();
    
    t2 = time(NULL);
    if(myRan == 0) printf("\nExecution took %f s\n", t2-t1);
#endif
    
    return 0;
}
