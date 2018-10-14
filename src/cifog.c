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
 
int cifog(confObj_t simParam, const double *redshift_list, grid_t *grid, sourcelist_t *sourcelist, const integral_table_t *integralTable, photIonlist_t *photIonBgList, const int num_cycles, const int myRank)
{
    double delta_redshift = 0.;
    int snap=-1;
    
    //-------------------------------------------------------------------------------
    // start the loop
    //-------------------------------------------------------------------------------
    
    for(int cycle=0; cycle<num_cycles; cycle++)
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
        
        cifog_step(simParam, grid, sourcelist, integralTable, photIonBgList, cycle, snap, myRank);
    }
    
    return EXIT_SUCCESS;
}

int cifog_step(confObj_t simParam, grid_t *grid, sourcelist_t *sourcelist, const integral_table_t *integralTable, photIonlist_t *photIonBgList, const int cycle, int snap, const int myRank)
{
    const double f = simParam->f;
    char photHIFile[MAXLENGTH], XionFile[MAXLENGTH];
    
    if(myRank==0)
    {
        printf("\n******************\nSNAP %d\t CYCLE %d\n******************\n", snap, cycle);
    }
    
    if(myRank==0) printf("\n++++\nreading sources/nion file for snap = %d... ", snap);
    read_update_nion(simParam, sourcelist, grid, snap);
    if(myRank==0) printf("done\n+++\n");
    
    if(simParam->solve_He == 1)
    {
        if(myRank==0) printf("\n++++\nreading sources/nion file for snap = %d... ", snap);
        read_update_nion_HeI(simParam, sourcelist, grid, snap);
        if(myRank==0) printf("done\n+++\n");
        
        if(myRank==0) printf("\n++++\nreading sources/nion file for snap = %d... ", snap);
        read_update_nion_HeII(simParam, sourcelist, grid, snap);
        if(myRank==0) printf("done\n+++\n");
    }
    
    if(myRank==0) printf("\n++++\nreading igm density file for snap = %d... ", snap);
    read_update_igm_density(simParam, grid, snap);
    if(myRank==0) printf("done\n+++\n");
  
    if(myRank==0) printf("\n++++\nreading igm clump file for snap = %d... ", snap);
    read_update_igm_clump(simParam, grid, snap);
    if(myRank==0) printf("done\n+++\n");
    
    //------------------------------------------------------------------------------
    // compute web model
    //------------------------------------------------------------------------------
    
    if(simParam->use_web_model == 1)
    {
        if(cycle == 0)
        {
            /* ----------------------------------------- */
            /* photoionization rate is homogeneous       */
            /* ----------------------------------------- */
            if(simParam->photHI_model == 0)
            {
                //set photoionization rate on grid to background value
                if(myRank==0) printf("\n++++\nsetting photoionization rate to background value... ");
                set_value_to_photoionization_field(grid, simParam);
                if(myRank==0) printf("\n photHI_bg = %e s^-1\n", simParam->photHI_bg);
                if(myRank==0) printf("done\n+++\n");
            }
            
            /* ----------------------------------------------------------------------------------------- */
            /* photoionization rate is given by mean mfp and exp(-r/mfp)/r^2 around sources distribution */
            /* ----------------------------------------------------------------------------------------- */
            else if(simParam->photHI_model == 11)
            {
                //set photoionization value according to the given list
                set_value_to_photHI_bg(simParam, get_photHI_from_redshift(photIonBgList, simParam->redshift));
                
                if(simParam->calc_mfp == 1)
                {
                    if(myRank==0) printf("\n++++\ncompute mean free path... ");
                    simParam->mfp = f*simParam->box_size/(simParam->h * (1.+simParam->redshift))/grid->nbins;
                    if(myRank==0) printf("\n mfp = %e Mpc at z = %e\n", simParam->mfp, simParam->redshift);
                    if(myRank==0) printf("done\n+++\n");
                }
                
                //compute spatial photoionization rate according to source distribution
                if(myRank==0) printf("\n++++\ncompute photoionization rate... ");
                compute_photHI(grid, simParam, 1);
                if(myRank==0) printf("done\n+++\n");
            }
            
            /* ----------------------------------------------------------------------------------------- */
            /* photoionization rate is given by mean mfp and exp(-r/mfp)/r^2 around sources distribution */
            /* ----------------------------------------------------------------------------------------- */
            else if(simParam->photHI_model == 1)
            {                
                if(simParam->calc_mfp == 1)
                {
                    if(myRank==0) printf("\n++++\ncompute mean free path... ");
                    simParam->mfp = f*simParam->box_size/(simParam->h * (1.+simParam->redshift))/grid->nbins;
                    if(myRank==0) printf("\n mfp = %e Mpc at z = %e\n", simParam->mfp, simParam->redshift);
                    if(myRank==0) printf("done\n+++\n");
                }
                
                //compute spatial photoionization rate according to source distribution
                if(myRank==0) printf("\n++++\ncompute photoionization rate... ");
                compute_photHI(grid, simParam, 0);
                if(myRank==0) printf("done\n+++\n");
            }
            
            /* ------------------------------------------------------- */
            /* photoionization rate is given by mfp of ionized regions */
            /* ------------------------------------------------------- */
            else if(simParam->photHI_model == 2)
            {                
                if(myRank==0) printf("\n++++\nset photoionization rate according to ionized regions... ");
                if(cycle != 0){
                    compute_photHI_ionizedRegions(grid, simParam);
                }else{
                    set_value_to_photoionization_field(grid,simParam);
                }
                if(myRank==0) printf("done\n+++\n");
            }
            
            else
            {
                if(myRank==0)
                {
                    printf("\n+++\nno supported photoionization rate model. Photoionization model is required for the web model. Abborting...\n");
                    exit(EXIT_FAILURE);
                }
            }
        }
        
        /* ------------------------------------------------------- */
        /* compute HI fraction (web model)                         */
        /* ------------------------------------------------------- */
        if(myRank==0) printf("\n++++\napply web model... ");
        compute_web_ionfraction(grid, simParam);
        if(myRank==0) printf("done\n+++\n");
        
        
        /* ------------------------------------------------------- */
        /* DISABLED: compute local mfp in each cell                */
        /* ------------------------------------------------------- */
        if(simParam->calc_mfp == -1)
        {
            //compute mean free paths
            if(myRank==0) printf("\n++++\ncompute mean free paths... ");
            compute_web_mfp(grid, simParam);
            if(myRank==0) printf("done\n+++\n");
        }
    }
    
    /* ---------------------------------------------------------------- */
    /* compute recombinations based on local photion rate & HI fraction */
    /* ---------------------------------------------------------------- */
    if(simParam->calc_recomb == 1 && simParam->const_recomb == 0)
    {
        //compute number of recombinations (HII, HeII & HeIII)
        if(myRank==0) printf("\n++++\ncompute number of recombinations... ");
        compute_number_recombinations(grid, simParam);
        if(myRank==0) printf("done\n+++\n");
    }
    else if(simParam->calc_recomb == 2 && simParam->const_recomb == 0)
    {
        //compute number of recombinations (HII, HeII & HeIII)
        if(myRank==0) printf("\n++++\ncompute number of recombinations... ");
        compute_number_recombinations_M2000(grid, simParam, simParam->recomb_table, integralTable);
        if(myRank==0) printf("done\n+++\n");
    }
    
    //------------------------------------------------------------------------------
    // compute number of recombinations for a constant rate
    //------------------------------------------------------------------------------
    
    if(simParam->const_recomb == 1)
    {
        //compute number of recombinations (HII, HeII & HeIII)
        if(myRank==0) printf("\n++++\ncompute number of recombinations... ");
        compute_number_recombinations_const(grid, simParam);
        if(myRank==0) printf("done\n+++\n");
    }

    //--------------------------------------------------------------------------------
    // apply tophat filter
    //--------------------------------------------------------------------------------
    
    //compute fraction Q
    if(myRank==0) printf("\n++++\nHII: computing relation between number of ionizing photons and absorptions... ");
    compute_cum_values(grid, simParam, 0);
    if(myRank==0) printf("done\n+++\n");
    
    //apply filtering
    if(myRank==0) printf("\n++++\nHII: apply tophat filter routine for ionization field... ");
    compute_ionization_field(simParam, grid, photIonBgList, 0);
    if(myRank==0) printf("done\n+++\n");
    
    //write ionization field to file
    for(int i=0; i<MAXLENGTH; i++) XionFile[i] = '\0';
    sprintf(XionFile, "%s_%02d", simParam->out_XHII_file, cycle);
    if(myRank==0) printf("\n++++\nwriting HI ionization field to file %s ... ", XionFile);
    save_to_file(grid->XHII, grid, XionFile);
    if(myRank==0) printf("done\n+++\n");
    
    if(simParam->solve_He == 1)
    {
        //compute fraction Q
        if(myRank==0) printf("\n++++\nHeII/HeIII: computing relation between number of ionizing photons and absorptions... ");
        compute_cum_values(grid, simParam, 1);
        compute_cum_values(grid, simParam, 2);
        if(myRank==0) printf("done\n+++\n");
    
        //apply filtering
        if(myRank==0) printf("\n++++\nHeII/HeIII: apply tophat filter routine for ionization field... ");
        compute_ionization_field(simParam, grid, photIonBgList, 1);
        compute_ionization_field(simParam, grid, photIonBgList, 2);
        if(myRank==0) printf("done\n+++\n");
        
        //write ionization field to file
        for(int i=0; i<MAXLENGTH; i++) XionFile[i] = '\0';
        sprintf(XionFile, "%s_%02d", simParam->out_XHeII_file, cycle);
        if(myRank==0) printf("\n++++\nwriting HeI ionization field to file %s ... ", XionFile);
        save_to_file(grid->XHeII, grid, XionFile);
        if(myRank==0) printf("done\n+++\n");
    
        //write ionization field to file
        for(int i=0; i<MAXLENGTH; i++) XionFile[i] = '\0';
        sprintf(XionFile, "%s_%02d", simParam->out_XHeIII_file, cycle);
        if(myRank==0) printf("\n++++\nwriting HeII ionization field to file %s ... ", XionFile);
        save_to_file(grid->XHeIII, grid, XionFile);
        if(myRank==0) printf("done\n+++\n");
    }
    
    if(simParam->use_web_model == 1 && simParam->write_photHI_file == 1)
    {
        //write photoionization rate field to file
        for(int i=0; i<MAXLENGTH; i++) photHIFile[i] = '\0';
        sprintf(photHIFile, "%s_%02d", simParam->out_photHI_file, cycle);
        if(myRank==0) printf("\n++++\nwriting HI photoionization field to file... ");
        save_to_file(grid->photHI, grid, photHIFile);
        if(myRank==0) printf("done\n+++\n");
    }
    
    return EXIT_SUCCESS;
}
