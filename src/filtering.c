#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "utils.h"
#include "confObj.h"
#include "grid.h"
#include "sources.h"

#include "fraction_q.h"
#include "convolution_fftw.h"

#include "self_shielding.h"
#include "density_distribution.h"

#define SQR(X) ((X) * (X))

void determine_mfp_nion(fftw_complex *frac_Q_smooth, fftw_complex *nion_smooth, fftw_complex *mfp_nion, int nbins, ptrdiff_t local_n0, double scale)
{
    double mfp_nion_tmp = 0.;
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                mfp_nion_tmp = creal(mfp_nion[i*nbins*nbins+j*nbins+k]);
                if(mfp_nion_tmp==0.)
                {
                    if(creal(frac_Q_smooth[i*nbins*nbins+j*nbins+k])>=1.)
                    {
                        mfp_nion[i*nbins*nbins+j*nbins+k] = scale*creal(nion_smooth[i*nbins*nbins+j*nbins+k])+0.*I;

                    }else{
                        mfp_nion[i*nbins*nbins+j*nbins+k] = creal(nion_smooth[i*nbins*nbins+j*nbins+k])/(double)nbins+0.*I;
                    }
                }else{
                    if(scale*creal(nion_smooth[i*nbins*nbins+j*nbins+k]) > mfp_nion_tmp && creal(frac_Q_smooth[i*nbins*nbins+j*nbins+k])>=1.)
                    {
                        mfp_nion[i*nbins*nbins+j*nbins+k] = scale*creal(nion_smooth[i*nbins*nbins+j*nbins+k])+0.*I;
                    }
                }
            }
        }
    }    
}

void determine_mfp(fftw_complex *frac_Q_smooth, fftw_complex *mfp, int nbins, ptrdiff_t local_n0, double scale)
{
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                if(creal(mfp[i*nbins*nbins+j*nbins+k])==0.){
                    if(creal(frac_Q_smooth[i*nbins*nbins+j*nbins+k])>=1.)
                    {
                        mfp[i*nbins*nbins+j*nbins+k] = scale + 0.*I;
                    }
                }
            }
        }
    }    
}

double determine_mean_mfp(fftw_complex *mfp, int nbins, ptrdiff_t local_n0)
{
    double sum_mfp=0.;
    int sum_cell=0;
    double result = 0.;
#ifdef __MPI
    double sum_mfp_all = 0.;
    int sum_cell_all = 0;
#endif
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                if(creal(mfp[i*nbins*nbins+j*nbins+k])>0.)
                {
                    sum_mfp += creal(mfp[i*nbins*nbins+j*nbins+k]);
                    sum_cell++;
                }
            }
        }
    }    
    
#ifdef __MPI
    MPI_Allreduce(&sum_mfp, &sum_mfp_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&sum_cell, &sum_cell_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    result = sum_mfp_all/((double)sum_cell_all);
    if(sum_cell_all <= 0) return 1./(double)nbins;
    else return result;
#else
    result = sum_mfp/((double)sum_cell);
    if(sum_cell <= 0) return 1./(double)nbins;
    else return result;
#endif
}

void construct_tophat_filter(fftw_complex *filter, int nbins, ptrdiff_t local_0_start, ptrdiff_t local_n0, float smooth_scale)
{
    double normCoeff;
    const double half_nbins = nbins*0.5;
    const double sq_smooth_scale = 0.25*SQR(smooth_scale);
    
    normCoeff = 0;
    for(int i=0; i<nbins; i++){      
        const double i_expr = half_nbins-fabs(i - half_nbins);        
        const double sq_i_expr = SQR(i_expr);
        for(int j=0; j<nbins; j++){
            const double j_expr = half_nbins - fabs(j - half_nbins);
            const double sq_j_expr = SQR(j_expr);
            for(int k=0; k<nbins; k++){
              const double k_expr = half_nbins - fabs(k - half_nbins);
              const double sq_k_expr = SQR(k_expr);
              
              const double expr = sq_i_expr + sq_j_expr + sq_k_expr - sq_smooth_scale;
              normCoeff += (expr<=0.) ? 1.:0.;
            }
        }
    }
    normCoeff = 1./normCoeff;
    
    for(int i=0; i<local_n0; i++)
    {
        const double i_expr = half_nbins - fabs(i + local_0_start - half_nbins);        
        const double sq_i_expr = SQR(i_expr);
        for(int j=0; j<nbins; j++)
        {
            const double j_expr = half_nbins - fabs(j - half_nbins);
            const double sq_j_expr = SQR(j_expr);
            for(int k=0; k<nbins; k++)
            {            
                const double k_expr = half_nbins - fabs(k - half_nbins);
                const double sq_k_expr = SQR(k_expr);
              
                const double expr = sq_i_expr + sq_j_expr + sq_k_expr - sq_smooth_scale;
                
                filter[i*nbins*nbins+j*nbins+k] = (expr<=0.) ? normCoeff+0.*I : 0.+0.*I;
            }
        }
    }
}

// versatile function (thisGrid is only used for grid dimensions and domain decomposition
void determine_ion_fractions(fftw_complex *cum_nion_smooth, int nbins, ptrdiff_t local_n0, int smallest_scale)
{
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                if(creal(cum_nion_smooth[i*nbins*nbins+j*nbins+k])>=1.)
                {
                    cum_nion_smooth[i*nbins*nbins+j*nbins+k] = 1.+0.*I;
                }else{
                    if(smallest_scale==0)
                    {
                        cum_nion_smooth[i*nbins*nbins+j*nbins+k] = 0.+0.*I;
                    }else{
                        if(creal(cum_nion_smooth[i*nbins*nbins+j*nbins+k])<1.e-8)
                        {
                            cum_nion_smooth[i*nbins*nbins+j*nbins+k] = 0.+0.*I;
                        }
                    }
                }
            }
        }
    }
}

void map_central_ionized_cell_to_sphere(fftw_complex *new_cum_nion_smooth, fftw_complex *cum_nion_smooth, fftw_complex *filter, grid_t *thisGrid)
{
    int nbins = thisGrid->nbins;
    int local_n0 = thisGrid->local_n0;
    double threshold = 1./((float)nbins*(float)nbins*(float)nbins);
    
    convolve_fft(thisGrid, filter, new_cum_nion_smooth, cum_nion_smooth);
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                if(creal(new_cum_nion_smooth[i*nbins*nbins+j*nbins+k])>threshold)
                {
                    new_cum_nion_smooth[i*nbins*nbins+j*nbins+k] = 1. + 0.*I;
                }else{
                    new_cum_nion_smooth[i*nbins*nbins+j*nbins+k] = 0. + 0.*I;
                }
            }
        }
    }
}

void map_central_ionized_cell_with_source_to_sphere(fftw_complex *new_cum_nion_smooth, fftw_complex *cum_nion_smooth, fftw_complex *filter, fftw_complex *cum_nion, grid_t *thisGrid)
{
    int nbins = thisGrid->nbins;
    int local_n0 = thisGrid->local_n0;
    double threshold = 1./((float)nbins*(float)nbins*(float)nbins);
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                if(creal(cum_nion[i*nbins*nbins+j*nbins+k]) <=0.)
                {
                    cum_nion_smooth[i*nbins*nbins+j*nbins+k] =  0. + 0.*I;
                }
            }
        }
    }
    
    convolve_fft(thisGrid, filter, new_cum_nion_smooth, cum_nion_smooth);
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                if(creal(new_cum_nion_smooth[i*nbins*nbins+j*nbins+k])>threshold)
                {
                    new_cum_nion_smooth[i*nbins*nbins+j*nbins+k] = 1. + 0.*I;
                }else{
                    new_cum_nion_smooth[i*nbins*nbins+j*nbins+k] = 0. + 0.*I;
                }
            }
        }
    }
}

// versatile function (thisGrid is only used for grid dimensions and domain decomposition
void choose_ion_fraction(fftw_complex *cum_nion_smooth, fftw_complex *Xion_tmp, grid_t *thisGrid)
{
    int nbins;
    ptrdiff_t local_n0;
    
    nbins = thisGrid->nbins;
    local_n0 = thisGrid->local_n0;
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                if(creal(cum_nion_smooth[i*nbins*nbins+j*nbins+k])>=creal(Xion_tmp[i*nbins*nbins+j*nbins+k]))
                {
                    Xion_tmp[i*nbins*nbins+j*nbins+k] =  cum_nion_smooth[i*nbins*nbins+j*nbins+k];
                }
            }
        }
    }
}

// versatile function (thisGrid is only used for grid dimensions and domain decomposition
void map_bubbles_to_nrec(fftw_complex *Xion_tmp, fftw_complex *nrec, grid_t *thisGrid)
{
      int nbins;
    ptrdiff_t local_n0;
    
    nbins = thisGrid->nbins;
    local_n0 = thisGrid->local_n0;
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                nrec[i*nbins*nbins+j*nbins+k] = nrec[i*nbins*nbins+j*nbins+k]*Xion_tmp[i*nbins*nbins+j*nbins+k];
            }
        }
    }
}

// versatile function (thisGrid is only used for grid dimensions and domain decomposition
void combine_bubble_and_web_model(fftw_complex *Xion_tmp, fftw_complex *Xion, grid_t *thisGrid)
{
    int nbins;
    ptrdiff_t local_n0;
    
    nbins = thisGrid->nbins;
    local_n0 = thisGrid->local_n0;
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                Xion[i*nbins*nbins+j*nbins+k] = Xion[i*nbins*nbins+j*nbins+k]*Xion_tmp[i*nbins*nbins+j*nbins+k];
            }
        }
    }
}

// versatile function (thisGrid is only used for grid dimensions and domain decomposition
void copy_grid_array(fftw_complex *Xion_tmp, fftw_complex *Xion, grid_t *thisGrid)
{
      int nbins;
    ptrdiff_t local_n0;
    
    nbins = thisGrid->nbins;
    local_n0 = thisGrid->local_n0;
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                Xion[i*nbins*nbins+j*nbins+k] = Xion_tmp[i*nbins*nbins+j*nbins+k];
            }
        }
    }
}

void update_web_model(grid_t *thisGrid, confObj_t simParam)
{
    const double f = simParam->f;

    thisGrid->mean_photHI = 0.;
    if(simParam->photHI_model == 2) compute_photHI_ionizedRegions(thisGrid, simParam);
    if(simParam->photHI_model == 1)
    {
        if(simParam->calc_mfp == 1)
        {
            set_mfp_Miralda2000(simParam);
            printf("\n M2000: mfp(photHI = %e) = %e Mpc at z = %e", simParam->photHI_bg, simParam->mfp, simParam->redshift);
            if(f*thisGrid->mean_mfp < simParam->mfp || simParam->photHI_bg < 1.e-12)
            {
                simParam->mfp = f*thisGrid->mean_mfp;
            }
            printf("\n mfp = %e Mpc at z = %e", simParam->mfp, simParam->redshift);
        }
        compute_photHI(thisGrid, simParam, 0);
    }
    compute_web_ionfraction(thisGrid, simParam);
}

void adapt_HeII_to_HeIII(grid_t *thisGrid)
{
    int nbins;
    ptrdiff_t local_n0;
    
    nbins = thisGrid->nbins;
    local_n0 = thisGrid->local_n0;
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                if(creal(thisGrid->XHeIII[i*nbins*nbins+j*nbins+k])==1.)
                {
                    thisGrid->XHeII[i*nbins*nbins+j*nbins+k] = 1.-creal(thisGrid->XHeIII[i*nbins*nbins+j*nbins+k]) + 0.*I;
                }
            }
        }
    }
}

void compute_ionization_field(confObj_t simParam, grid_t *thisGrid, int specie)
{
    int nbins;
    float smooth_scale;
    float box_size;
    float lin_scales, inc_log_scales, max_scale;
    
#ifdef DEBUG_FILTERING
    char scale_string[8];
    char Q_file[MAXLENGTH];
#endif
    
    fftw_complex *filter;
    fftw_complex *cum_nion_smooth;
    fftw_complex *cum_nabs_smooth;
    fftw_complex *frac_Q_smooth;
    fftw_complex *Xion_tmp;
    fftw_complex *nion_smooth = NULL;
    fftw_complex *mfp_tmp = NULL;    
    
    fftw_complex *cum_nion;
    fftw_complex *cum_nabs;
    fftw_complex *cum_nrec;
    fftw_complex *Xion;
    fftw_complex *nion = NULL;
        
    if(specie == 1){
        cum_nion = thisGrid->cum_nion_HeI;
        cum_nabs = thisGrid->cum_nabs_HeI;
        cum_nrec = thisGrid->cum_nrec_HeI;
        Xion = thisGrid->XHeII;
    }else if(specie == 2){
        cum_nion = thisGrid->cum_nion_HeII;
        cum_nabs = thisGrid->cum_nabs_HeII;
        cum_nrec = thisGrid->cum_nrec_HeII;
        Xion = thisGrid->XHeIII;
    }else{
        cum_nion = thisGrid->cum_nion;
        cum_nabs = thisGrid->cum_nabs;
        cum_nrec = thisGrid->cum_nrec;
        Xion = thisGrid->XHII;
        nion = thisGrid->nion;  //only needed to compute mfp
    }
    
#ifdef __MPI
    ptrdiff_t alloc_local, local_n0, local_0_start;
#else
    ptrdiff_t local_n0, local_0_start;
#endif
    
    float inc, factor_exponent, lin_bins;
    int num_scales;

    nbins = thisGrid->nbins;
    box_size = thisGrid->box_size;
    lin_scales = thisGrid->lin_scales;
    inc_log_scales = thisGrid->inc_log_scales;
    max_scale = thisGrid->max_scale;
    
#ifdef __MPI
    alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
    filter = fftw_alloc_complex(alloc_local);
    cum_nion_smooth = fftw_alloc_complex(alloc_local);
    cum_nabs_smooth = fftw_alloc_complex(alloc_local);
    frac_Q_smooth = fftw_alloc_complex(alloc_local);
    
    Xion_tmp = fftw_alloc_complex(alloc_local);
    mfp_tmp = fftw_alloc_complex(alloc_local);
    if(simParam->photHI_model == 2)
    {
        nion_smooth = fftw_alloc_complex(alloc_local);
    }
#else
    local_0_start = thisGrid->local_0_start;
    local_n0 = thisGrid->local_n0;
    filter = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    if(filter == NULL)
    {
        fprintf(stderr, "filter in compute_ionization_field (filtering.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    cum_nion_smooth = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    if(cum_nion_smooth == NULL)
    {
        fprintf(stderr, "cum_nion_smooth in compute_ionization_field (filtering.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    cum_nabs_smooth = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    if(cum_nabs_smooth == NULL)
    {
        fprintf(stderr, "cum_nabs_smooth in compute_ionization_field (filtering.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    frac_Q_smooth = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    if(frac_Q_smooth == NULL)
    {
        fprintf(stderr, "frac_Q_smooth in compute_ionization_field (filtering.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    Xion_tmp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    if(Xion_tmp == NULL)
    {
        fprintf(stderr, "Xion_tmp in compute_ionization_field (filtering.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    mfp_tmp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    if(mfp_tmp == NULL)
    {
        fprintf(stderr, "mfp_tmp in compute_ionization_field (filtering.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    if(simParam->photHI_model == 2) 
    {
        nion_smooth = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        if(nion_smooth == NULL)
        {
            fprintf(stderr, "nion_smooth in compute_ionization_field (filtering.c) could not be allocated\n");
            exit(EXIT_FAILURE);
        }
    }
#endif

    initialize_grid(Xion_tmp, nbins, local_n0, 0.);
    initialize_grid(mfp_tmp, nbins, local_n0, 0.);
    if(simParam->photHI_model == 2 && (specie != 1 && specie !=2)) 
    {
        initialize_grid(nion_smooth, nbins, local_n0, 0.);
    }
    
    /* ---------------------------------------------- */
    /* compute number and sizes of smoothing scales   */
    /* ---------------------------------------------- */
    lin_bins = lin_scales/box_size*(float)nbins;
    factor_exponent = inc_log_scales/lin_bins;
    num_scales = (log(max_scale/lin_scales) + lin_bins*log(1.+factor_exponent))/log(1. + factor_exponent );
    printf("\n #linear bins = %e\t factor_exponent = %e\t #tophat filter sizes = %d\n", lin_bins, factor_exponent, num_scales);
    
    
    /* ----------------------------------------- */
    /* loop over different smoothing scales      */
    /* ----------------------------------------- */
    for(int scale=0; scale<num_scales; scale++)
    {
        /* ----------------------------------------- */
        /* map Nion and Nabs value to temporal grids */
        /* ----------------------------------------- */
        for(int i=0; i<local_n0; i++)
        {
            for(int j=0; j<nbins; j++)
            {
                for(int k=0; k<nbins; k++)
                {
                    cum_nion_smooth[i*nbins*nbins+j*nbins+k] = creal(cum_nion[i*nbins*nbins+j*nbins+k])+0.*I;
                    cum_nabs_smooth[i*nbins*nbins+j*nbins+k] = creal(cum_nabs[i*nbins*nbins+j*nbins+k])+0.*I;
                    if(simParam->photHI_model == 2)
                    {
                        nion_smooth[i*nbins*nbins+j*nbins+k] = creal(nion[i*nbins*nbins+j*nbins+k])+0.*I;
                    }
                }
            }
        }
        
        /* ----------------------------------------- */
        /* compute smoothing scale                   */
        /* ----------------------------------------- */
        inc = (float)num_scales - (float)scale;
        
        if(inc <= lin_bins) smooth_scale = inc;
        else smooth_scale = lin_bins*pow(1. + factor_exponent, inc-lin_bins);
        
        printf("  inc = %e\t scale = %d\t smooth_scale = %e bins or %e cMpc\n",inc, scale, smooth_scale, (float)smooth_scale/(float)nbins*box_size);

        /* -------------------------------------------- */
        /* coonstruct tophat filter for smoothing scale */
        /* -------------------------------------------- */
        construct_tophat_filter(filter, nbins, local_0_start, local_n0, smooth_scale);
        
        /* -------------------------------------------- */
        /* convolution of tophat filter and Nion & Nabs */
        /* -------------------------------------------- */
        convolve_fft(thisGrid, filter, cum_nion_smooth, cum_nion);
        convolve_fft(thisGrid, filter, cum_nabs_smooth, cum_nabs);
        if(simParam->photHI_model == 2) convolve_fft(thisGrid, filter, nion_smooth, nion);

        /* -------------------------------------------------------- */
        /* compute fraction of number of ionization and absorptions */
        /* -------------------------------------------------------- */ 
        compute_Q(thisGrid, frac_Q_smooth, cum_nion_smooth, cum_nabs_smooth);
//         save_to_file(frac_Q_smooth, thisGrid, "Q.dat");
        
        /* -------------------------------------------- */
        /* derive ionized regions                       */
        /* -------------------------------------------- */ 
        if(scale==num_scales-1)
        {
            determine_ion_fractions(frac_Q_smooth, nbins, local_n0, 1);
        }else{
            determine_ion_fractions(frac_Q_smooth, nbins, local_n0, 0);
        }
        
#ifdef DEBUG_FILTERING
        sprintf(scale_string, "%03d", scale);
        for(int i=0; i<MAXLENGTH; i++) Q_file[i] = '\0';
        strcat(Q_file, "Q1_");
        strcat(Q_file, scale_string);
        strcat(Q_file, ".dat");
        save_to_file(frac_Q_smooth, thisGrid, Q_file);
#endif
        
        if(simParam->ionize_sphere == 1 && scale != num_scales-1)
        {
            map_central_ionized_cell_to_sphere(frac_Q_smooth, frac_Q_smooth, filter, thisGrid);
        }
        else if(simParam->ionize_sphere == 2 && scale != num_scales-1)
        {
            map_central_ionized_cell_with_source_to_sphere(frac_Q_smooth, frac_Q_smooth, filter, cum_nion, thisGrid);
        }
        
#ifdef DEBUG_FILTERING
        sprintf(scale_string, "%03d", scale);
        for(int i=0; i<MAXLENGTH; i++) Q_file[i] = '\0';
        strcat(Q_file, "Q2_");
        strcat(Q_file, scale_string);
        strcat(Q_file, ".dat");
        save_to_file(frac_Q_smooth, thisGrid, Q_file);
#endif
        
        choose_ion_fraction(frac_Q_smooth, Xion_tmp, thisGrid);

        /* -------------------------------------------- */
        /* derive mfp from ionized regions              */
        /* -------------------------------------------- */ 
        if(specie == 0)
        {
            if(simParam->photHI_model == 2)
            {
                determine_mfp_nion(frac_Q_smooth, nion_smooth, mfp_tmp, nbins, local_n0, (double)smooth_scale/(double)nbins);
            }else{
                determine_mfp(frac_Q_smooth, mfp_tmp, nbins, local_n0, (double)smooth_scale/(double)nbins);
            }
        }
    }
    
    
    /* -------------------------------------------------------- */
    /* save mean mfp of ionization field                        */
    /* -------------------------------------------------------- */ 
    double phys_boxsize = simParam->box_size/(simParam->h*(1.+simParam->redshift));
    if(specie == 0)
    {
        if(simParam->photHI_model != 2)
        {
            thisGrid->mean_mfp = determine_mean_mfp(mfp_tmp, nbins, local_n0)*phys_boxsize;
        }else{
            thisGrid->mean_mfp = 0.;
        }
    }
    
    /* -------------------------------------------------------- */
    /* transfer mfp calculation to photHI grid                  */
    /* -------------------------------------------------------- */ 
    if(simParam->photHI_model == 2 && specie == 0)
    {
        copy_grid_array(mfp_tmp, thisGrid->photHI, thisGrid);
    }
    
    /* -------------------------------------------------------- */
    /* apply web model                                          */
    /* -------------------------------------------------------- */ 
    if(simParam->use_web_model == 1 && specie == 0)
    {        
        update_web_model(thisGrid, simParam);
        combine_bubble_and_web_model(Xion_tmp, Xion, thisGrid);
    }else{
        copy_grid_array(Xion_tmp, Xion, thisGrid);
    }
    map_bubbles_to_nrec(Xion_tmp, cum_nrec, thisGrid);
    
    /* -------------------------------------------------------- */
    /* combine XHeII and XHeIII                                 */
    /* -------------------------------------------------------- */ 
    if(specie == 2)
    {
        adapt_HeII_to_HeIII(thisGrid);
    }
    
    /* -------------------------------------------------------- */
    /* deallocation of arrays                                   */
    /* -------------------------------------------------------- */
    fftw_free(filter);
    fftw_free(cum_nion_smooth);
    fftw_free(cum_nabs_smooth);
    fftw_free(frac_Q_smooth);
    fftw_free(Xion_tmp);
    fftw_free(mfp_tmp);
    if(simParam->photHI_model == 2)
    {
        fftw_free(nion_smooth);
    }
}
