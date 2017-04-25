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

#include "cifog.h"

int main (int argc, /*const*/ char * argv[]) { 
#ifdef __MPI
    int size = 1;
#endif
    int myRank = 0;

    char iniFile[MAXLENGTH];
    confObj_t simParam;
    
    double *redshift_list = NULL;
    
    grid_t *grid = NULL;
    
    sourcelist_t *sourcelist = NULL;
    
    integral_table_t *integralTable = NULL;
    
    photIonlist_t *photIonBgList = NULL;
    
    double t1, t2;
    
    double zstart = 0., zend = 0., delta_redshift = 0.;
    int num_cycles;
        
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
    
    if(myRank==0)
    {
        printf("\n++++\nTEST OUTPUT\n");
    //     printf("densSS = %e\n", ss_calc_densSS(simParam, 1.e-13, 1.e4, 6.));
    //     printf("densSS = %e\n", ss_calc_densSS(simParam, 5.1e-11, 1.e4, 14.75));
    //     printf("densSS = %e\n", ss_calc_densSS(simParam, 1.e-12, 1.e4, 9.));
    //     printf("densSS = %e\n", ss_calc_densSS(simParam, 1.e-12, 1.e4, 7.)*simParam->omega_b*simParam->h*simParam->h*rho_g_cm/mp_g*8.*8.*8./(1.-simParam->Y));
        printf(" mean free paths for T=10^4K and photHI_bg = 5.e-13 :\n");
        printf(" z = 6: mfp = %e\n", calc_local_mfp(simParam, 1., 0.5e-12, 1.e4, 6.));
        printf(" z = 6: mfp(M2000) = %e\n", dd_calc_mfp(simParam, 0.5e-12, 1.e4, 6.));
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
        printf("z = 14.75: mfp(M2000) = %e\n", dd_calc_mfp(simParam, 0.5e-13, 1.e4, 14.75));
        printf("z = 14.75: mfp(M2000) = %e\n", dd_calc_mfp(simParam, 0.5e-14, 1.e4, 14.75));
        printf("z = 14.75: mfp(M2000) = %e\n", dd_calc_mfp(simParam, 0.5e-15, 1.e4, 14.75));
        printf("done\n+++\n");
    }
    
    //verify that helium runs contain helium!
    if(simParam->solve_He == 1)
    {
        if(simParam->Y <= 0.)
        {
            fprintf(stderr, "If you include helium its mass fraction Y should be larger than zero!\n");
            exit(EXIT_FAILURE);
        }
    }
    
    //read redshift files with outputs
    redshift_list = NULL;
    if(myRank==0) printf("\n++++\nreading redshift list of files and outputs... ");
    redshift_list = read_redshift_list(simParam->redshift_file, num_cycles);
    if(redshift_list != NULL) num_cycles = num_cycles - 1;
    if(myRank==0) printf("done\n+++\n");
    
    //read files (allocate grid)
    grid = initGrid();
    if(myRank==0) printf("\n++++\nreading files to grid... ");
    read_files_to_grid(grid, simParam);
    if(myRank==0) printf("done\n+++\n");
    
    //read photoionization background values 
    if(myRank==0) printf("\n++++\nreading photoionization background rates... ");
    photIonBgList = read_photIonlist(simParam->photHI_bg_file);
    if(myRank==0) printf("done\n+++\n");
    
    if(simParam->calc_recomb == 2)
    {
        //read table for recombinations
        if(myRank==0) printf("\n++++\nread table for recombinations... ");
        integralTable = initIntegralTable(simParam->zmin, simParam->zmax, simParam->dz, simParam->fmin, simParam->fmax, simParam->df, simParam->dcellmin, simParam->dcellmax, simParam->ddcell);
        if(myRank==0) printf("done\n+++\n");
    }

    printf("\nThis run computes %d times the ionization field (num_cycles)\n", num_cycles);
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
    
    cifog(simParam, redshift_list, grid, sourcelist, integralTable, photIonBgList, num_cycles, myRank);
    

    //--------------------------------------------------------------------------------
    // deallocating grids
    //--------------------------------------------------------------------------------
    if(simParam->calc_recomb == 2)
    {
        //read table for recombinations
        if(myRank==0) printf("\n++++\ndeallocating table for recominsations... ");
        free(integralTable);
        if(myRank==0) printf("done\n+++\n");
    }

    if(myRank==0) printf("\n++++\ndeallocating background photionization rate list... ");
    deallocate_photIonlist(photIonBgList);
    if(myRank==0) printf("done\n+++\n");

    //deallocate grid
    if(myRank==0) printf("\n++++\ndeallocating grid ...");
    deallocate_grid(grid, simParam);
    if(myRank==0) printf("done\n+++\n");
    
    //deallocate redshift list
    if(myRank==0) printf("\n++++\ndeallocating redshift list ...");
    deallocateRedshift_list(redshift_list);
    if(myRank==0) printf("done\n+++\n");
    
    confObj_del(&simParam);
    
    if(myRank==0) printf("\nFinished\n");
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
