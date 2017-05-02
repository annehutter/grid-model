#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>    //included because M_PI is not defined in <math.h>
#include <assert.h>
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

#include "convolution_fftw.h"
#include "redshift_tools.h"
#include "self_shielding.h"

#define SQR(X) ((X) * (X))
#define CUB(X) ((X) * (X) * (X))

double ss_calc_modPhotHI(double densH, double densSS)
{
    double quot;
    
    quot = densH/densSS;
    return 0.98*pow(1.+pow(quot,1.64),-2.28)+0.02*pow(1.+quot,-0.84);
}

double ss_calc_densSS(confObj_t simParam, double photHI, double temperature, double redshift)
{
    double tmp, tmp2, rec_rate;
    const double omega_b = simParam->omega_b;
    const double omega_m = simParam->omega_m;
    const double h = simParam->h;
    const double Y = simParam->Y;
    const double fg = omega_b/omega_m;
    const double mu = 4./(8.-5.*Y);
    double rho;
    
    if(simParam->default_mean_density == 1){
        rho = 3.*SQR(H0)/(8.*M_PI*G)/SQR(h);
    }else{
        rho = simParam->mean_density*mp_g;
    }
    
    tmp = mu*G/(M_PI*gamma_gas*boltzman_cgs);
    tmp2 = SQR(mp_g)/(fg*SQR(1.-Y)*SQR(1.-0.75*Y)*SQR(sigma_HI));
    rec_rate = recomb_HII*pow(temperature*1.e-4,-0.8);    //definition of the HII recombination rate (Fukugita & Kawasaki 1994)
    
    return pow(tmp*tmp2*SQR(photHI)/temperature/SQR(rec_rate),1./3.)/(omega_b*SQR(h)*rho/mp_g*CUB(1.+redshift));    //checked against Mesinger 2015!
    
}

double ss_calc_XHII(double dens, double photHI, double temp, double Y)
{
    double tmp, tmp2;
    
    tmp = photHI*(1.-Y)/((1.-0.75*Y)*dens*recomb_HII);
    tmp2 = 0.5*(tmp + 2. - sqrt((tmp +2.)*(tmp + 2.)-4.));
    if((1.-tmp2)>1.) return 1.;
    else return (1.-tmp2);
}

// double calc_photHI_source(source_t *thisSource, double mfp_inv, double boxsize_Mpc, float x, float y, float z)
// {
//     double r;
//     
//     const double dx = thisSource->pos[0]-x;
//     const double dy = thisSource->pos[1]-y;
//     const double dz = thisSource->pos[2]-z;
//     
//     r = sqrt(dx*dx+dy*dy+dz*dz)*boxsize_Mpc;
//     
//     return sigma_HI*thisSource->Nion*exp(-r*mfp_inv)/(r*r);
// }
// 
// void compute_photoionization_field(grid_t *thisGrid, sourcelist_t *thisSourcelist, confObj_t simParam)
// {
//     int numSources;
//     
//       int nbins;
//     double nbins_inv, boxsize_cm_inv;
//     ptrdiff_t local_0_start, local_n0;
//     
//     double mfp_inv;
//     
//     double tmp, sum;
//     
//     numSources = thisSourcelist->numSources;
//     
//     nbins = thisGrid->nbins;
//     nbins_inv = 1./(double)nbins;
//     boxsize_cm_inv = (1.+simParam->redshift)/(simParam->box_size*Mpc_cm);    //inverse boxsize at z in cm
// 
//     local_0_start = thisGrid->local_0_start;
//     local_n0 = thisGrid->local_n0;
//     
//     mfp_inv = 1./(simParam->mfp);
//     
//     sum = 0.;
//     for(int comz=0; comz<local_n0; comz++)
//     {
//         printf("z = %d\n",comz);
//           const float z = (comz+local_0_start)*nbins_inv;
//         for(int comy=0; comy<nbins; comy++)
//         {
//               const float y = comy*nbins_inv;
//             for(int comx=0; comx<nbins; comx++)
//             {
//                 const float x = comx*nbins_inv;
//                 tmp = 0.;
//                 const source_t *sources = thisSourcelist->source;
//                 for(int i=0; i<numSources; i++)
//                 {
//                     const source_t *thisSource = &sources[i];
//                     float dx = thisSource->pos[0]-x;
//                     float dy = thisSource->pos[1]-y;
//                     float dz = thisSource->pos[2]-z;
//                     
//                     if(dx>0.5) dx = dx - 1.0f;
//                     if(dy>0.5) dy = dy - 1.0f;
//                     if(dz>0.5) dz = dz - 1.0f;
//                     if(dx<-0.5) dx = dx + 1.0f;
//                     if(dy<-0.5) dy = dy + 1.0f;
//                     if(dz<-0.5) dz = dz + 1.0f;
//                     
//                     const float r2 = (dx*dx+dy*dy+dz*dz+0.00001f);
//                     const float r = sqrtf(r2);
//                     const float r2_inv = 1.0f/r2;
//                     
//                     tmp += sigma_HI*boxsize_cm_inv*thisSource->Nion*thisSource->fesc*boxsize_cm_inv*exp(-r*mfp_inv)*r2_inv;
// //                     printf("source: %e\t%e\t%e\n",tmp,r,exp(-r*mfp_inv));
// //                   tmp += calc_photHI_source(thisSourcelist->source[source], mfp_inv, boxsize_Mpc, x, y, z);
//                 }
//                 thisGrid->photHI[comz*nbins*nbins + comy*nbins + comx] = tmp + 0.*I;
//                 sum += tmp;
//             }
//         }
//     }
//     thisGrid->mean_photHI = sum;
// #ifdef __MPI
//     MPI_Allreduce(&sum, &thisGrid->mean_photHI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
// #endif
//     thisGrid->mean_photHI = thisGrid->mean_photHI*nbins_inv*nbins_inv*nbins_inv;
// }

void construct_photHI_filter(fftw_complex *filter, grid_t *thisGrid, confObj_t simParam)
{
    ptrdiff_t local_n0, local_0_start;
    int nbins;
    
    nbins = thisGrid->nbins;
    local_n0 = thisGrid->local_n0;
    local_0_start = thisGrid->local_0_start;
    
    const double half_nbins = nbins*0.5;
    const double mfp_inv = 1./simParam->mfp;   // in physical Mpc^-1
    const double factor = (thisGrid->box_size/simParam->h/(1.+simParam->redshift))/nbins;
    const double sq_factor = factor*factor;
    
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
              
                double expr = (sq_i_expr + sq_j_expr + sq_k_expr + 0.25)*sq_factor;
                
                filter[i*nbins*nbins+j*nbins+k] = exp(-sqrt(expr)*mfp_inv)/expr + 0.*I;
                
                if(creal(filter[i*nbins*nbins+j*nbins+k])<=0.) printf("%d: %e\n",(int)local_0_start,creal(filter[i*nbins*nbins+j*nbins+k]));
            }
        }
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void convolve_fft_photHI(grid_t *thisGrid, fftw_complex *filter, fftw_complex *nion_smooth, const double factor_photHI)
{
    int nbins;
    double factor;
    
#ifdef __MPI
    ptrdiff_t alloc_local, local_n0, local_0_start;
    ptrdiff_t n0, n1, n2, local_n1, local_1_start, local_n0x, local_0x_start;
#else
    ptrdiff_t local_n0;
#endif
    
    double mean_photHI;
    
    fftw_complex *nion;
    fftw_complex *nion_ft, *filter_ft;
    fftw_plan plan_nion, plan_filter, plan_back;
    
    nbins = thisGrid->nbins;
    nion = thisGrid->nion;
    local_n0 = thisGrid->local_n0;
    
#ifdef __MPI
    local_0_start = thisGrid->local_0_start;
    
    n0 = nbins;
    n1 = nbins;
    n2 = nbins;
    
    alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
    assert(local_0_start == thisGrid->local_0_start);
    assert(local_n0 == thisGrid->local_n0);
    
    nion_ft = fftw_alloc_complex(alloc_local);
    filter_ft = fftw_alloc_complex(alloc_local);
    
    plan_nion = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, nion, nion_ft, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_filter = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, filter, filter_ft, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);
#else 
    nion_ft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    filter_ft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    
    plan_nion = fftw_plan_dft_3d(nbins, nbins, nbins, nion, nion_ft, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_filter = fftw_plan_dft_3d(nbins, nbins, nbins, filter, filter_ft, FFTW_FORWARD, FFTW_ESTIMATE);
#endif
    
    fftw_execute(plan_nion);
    fftw_execute(plan_filter);

#ifdef __MPI
    fftw_mpi_local_size_3d_transposed(n0, n1, n2, MPI_COMM_WORLD, &local_n0x, &local_0x_start, &local_n1, &local_1_start);
#endif
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                nion_ft[i*nbins*nbins+j*nbins+k] = nion_ft[i*nbins*nbins+j*nbins+k]*filter_ft[i*nbins*nbins+j*nbins+k];
            }
        }
    }
    
#ifdef __MPI
    plan_back = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, nion_ft, nion_smooth, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);
#else
    plan_back = fftw_plan_dft_3d(nbins, nbins, nbins, nion_ft, nion_smooth, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
    fftw_execute(plan_back);
    
    factor = 1./(nbins*nbins*nbins);
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                nion_smooth[i*nbins*nbins+j*nbins+k] = factor*nion_smooth[i*nbins*nbins+j*nbins+k];
                if(creal(nion_smooth[i*nbins*nbins+j*nbins+k]) < 0.) nion_smooth[i*nbins*nbins+j*nbins+k] = 0. + 0.*I;
            }
        }
    }
    
    mean_photHI = creal(nion_smooth[0])*factor_photHI;
    
    fftw_destroy_plan(plan_nion);
    fftw_destroy_plan(plan_filter);
    fftw_destroy_plan(plan_back);
    
    fftw_free(nion_ft);
    fftw_free(filter_ft);
    
#ifdef __MPI
    MPI_Bcast(&mean_photHI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif    
    thisGrid->mean_photHI = mean_photHI;
}

void replace_convolve_fft_photHI(grid_t *thisGrid, confObj_t simParam, fftw_complex *nion_smooth, const double factor_photHI)
{
    int nbins = thisGrid->nbins;
    ptrdiff_t local_n0 = thisGrid->local_n0;
    
    double mean_photHI;
    
    const double mfp_inv = 1./simParam->mfp;
    const double factor = (thisGrid->box_size/simParam->h/(1.+simParam->redshift))/nbins;
    const double sq_factor = factor*factor;
    const double expr = 0.25 * sq_factor;
    const double factor_nion = exp(-sqrt(expr)*mfp_inv)/expr;
    
    double sum = 0.;
    int sum_int = 0;
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                nion_smooth[i*nbins*nbins+j*nbins+k] = creal(nion_smooth[i*nbins*nbins+j*nbins+k])*factor_nion + 0.*I;
                sum += creal(nion_smooth[i*nbins*nbins+j*nbins+k]);
                if(creal(nion_smooth[i*nbins*nbins+j*nbins+k]) > 0.)
                {
                    sum_int++;
//                     printf("nion = %e\n", creal(nion_smooth[i*nbins*nbins+j*nbins+k]));
                }
            }
        }
    }
    
    double mean_sep_cells = nbins/pow(sum_int,1./3.);
    const double expr_mean = mean_sep_cells * sq_factor;
    const double factor_mean = exp(-sqrt(expr_mean)*mfp_inv)/expr_mean;
    double value;
    if(sum_int > 0)
    {
        value = sum / (factor_nion * (double)sum_int);
    }else{
        value = 0.;
    }
        
    sum = 0.;
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                if(creal(nion_smooth[i*nbins*nbins+j*nbins+k]) <= 0.)
                {
                    nion_smooth[i*nbins*nbins+j*nbins+k] = factor_mean*value + 0.*I;
                }
                sum += creal(nion_smooth[i*nbins*nbins+j*nbins+k]);
            }
        }
    }
    
#ifdef __MPI
    MPI_Allreduce(&sum, &mean_photHI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    mean_photHI = sum;
#endif
    mean_photHI = mean_photHI/(nbins*nbins*nbins);
    mean_photHI = mean_photHI*factor_photHI;

    thisGrid->mean_photHI = mean_photHI;
}

void compute_photHI(grid_t *thisGrid, confObj_t simParam, int rescale)
{
#ifdef __MPI
    ptrdiff_t alloc_local, local_n0, local_0_start;
#else
    ptrdiff_t local_n0;
#endif
    int nbins = thisGrid->nbins;
    fftw_complex *filter;
    fftw_complex *nion;

    const double z = simParam->redshift;
    const double cellsize_phys = (thisGrid->box_size/simParam->h/(1.+z))/nbins;
    const double mfp_index = simParam->mfp/cellsize_phys;

    const double alpha = simParam->source_slope_index;
    const double beta = 3.;
    const double f = simParam->factor;
    const double factor_k = 4.*M_PI*CUB(f)*0.632121;
    const double factor_photHI = sigma_HI*alpha/((alpha+beta)*SQR(Mpc_cm)*factor_k);
            
    if(thisGrid->mean_photHI == 0.)
    {
#ifdef __MPI
        alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
        filter = fftw_alloc_complex(alloc_local);
        nion = fftw_alloc_complex(alloc_local);
        assert(local_n0 == thisGrid->local_n0);
        assert(local_0_start == thisGrid->local_0_start);
#else
        local_n0 = thisGrid->local_n0;
        filter = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
        nion = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
#endif
                
        if(simParam->mfp > cellsize_phys)
        {
            //construct exp(-r/mfp)/(r*r) filter
            construct_photHI_filter(filter, thisGrid, simParam);
        }
                
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        
        for(int i=0; i<local_n0; i++)
        {
            for(int j=0; j<nbins; j++)
            {
                for(int k=0; k<nbins; k++)
                {
                    nion[i*nbins*nbins+j*nbins+k] = creal(thisGrid->nion[i*nbins*nbins+j*nbins+k]) + 0.*I;
                }
            }
        }
        
        if(simParam->mfp > cellsize_phys)
        {
            //apply filter to Nion field
            convolve_fft_photHI(thisGrid, filter, nion, factor_photHI);
        }
        else
        {
            printf("\n mfp too small to do fft mapping...\t");
            replace_convolve_fft_photHI(thisGrid, simParam, nion, factor_photHI);
        }

        double rescale_factor = 1.;
        if(rescale == 1)
        {
            rescale_factor = simParam->photHI_bg/thisGrid->mean_photHI;
            printf("\n fitting the photoionization rate field with mean value %e to the given background value by multiplying with a factor = %e\n", thisGrid->mean_photHI, rescale_factor);
        }
                
        for(int i=0; i<local_n0; i++)
        {
            for(int j=0; j<nbins; j++)
            {
                for(int k=0; k<nbins; k++)
                {
                    if(creal(nion[i*nbins*nbins+j*nbins+k]) < 0.) printf("\n photHI = %e < 0. !!!", creal(nion[i*nbins*nbins+j*nbins+k]));
                    thisGrid->photHI[i*nbins*nbins+j*nbins+k] = creal(nion[i*nbins*nbins+j*nbins+k])*factor_photHI*rescale_factor + 0.*I;
                }
            }
        }

        if(simParam->padded_box != 0.)
        {
            thisGrid->mean_photHI = thisGrid->mean_photHI*simParam->padded_box;
        }
        simParam->photHI_bg = thisGrid->mean_photHI;
                
        fftw_free(filter);
        fftw_free(nion);
    }

    printf("\n mfp = %e Mpc\tmfp_index = %e cells\tfactork = %e", simParam->mfp, mfp_index, factor_k);
    printf("\n mean photHI (accounting for all cells) = %e", simParam->photHI_bg);
    printf("\n actual mean photHI (accounting only for ionized regions) = %e\n",calc_mean_photoionization_ionized_field(thisGrid));
}

void set_value_to_photoionization_field(grid_t *thisGrid, confObj_t simParam)
{
    ptrdiff_t local_n0;
    int nbins;
    
    local_n0 = thisGrid->local_n0;
    nbins = thisGrid->nbins;
    
    initialize_grid(thisGrid->photHI, nbins, local_n0, simParam->photHI_bg);
    thisGrid->mean_photHI = simParam->photHI_bg;
}

void set_value_to_photHI_bg(grid_t *thisGrid, confObj_t simParam, double value)
{
    thisGrid->mean_photHI = value;
    simParam->photHI_bg = value;
}

void compute_photHI_ionizedRegions(grid_t *thisGrid, confObj_t simParam)
{
    ptrdiff_t local_n0 = thisGrid->local_n0;
    int nbins = thisGrid->nbins;
    
    double z = simParam->redshift;
    double alpha = simParam->source_slope_index;
    double beta = 3.;
    
    double factor = thisGrid->box_size/simParam->h/(1.+z)*Mpc_cm;
    double len_cell = factor/nbins;         //physical length of a cell in cm
    double factor2 = 2.*sigma_HI*alpha/((alpha+beta)*len_cell*len_cell*len_cell); 
    
    if(thisGrid->mean_photHI == 0.)
    {
        for(int i=0; i<local_n0; i++)
        {
            for(int j=0; j<nbins; j++)
            {
                for(int k=0; k<nbins; k++)
                {
                    thisGrid->photHI[i*nbins*nbins+j*nbins+k] = creal(thisGrid->photHI[i*nbins*nbins+j*nbins+k])*factor*factor2 + 0.*I;
                }
            }
        }
        thisGrid->mean_photHI = get_mean_grid(thisGrid->photHI, nbins, local_n0);
    }
}

void compute_web_ionfraction(grid_t *thisGrid, confObj_t simParam)
{
      int nbins;
    ptrdiff_t local_n0;
    int cell;
    
    double mean_numdensity_H;
    double photHI;
    double densSS;
    double mod_photHI;
    double redshift;
    double temperature;
    double correct_HeII;
    
    nbins = thisGrid->nbins;
    local_n0 = thisGrid->local_n0;
    
    redshift = simParam->redshift;
    temperature = 1.e4;
    mean_numdensity_H = simParam->mean_density*(1.+redshift)*(1.+redshift)*(1.+redshift);

    if(simParam->default_mean_density == 1){
        mean_numdensity_H = rho_g_cm/mp_g*simParam->h*simParam->h*simParam->omega_b*(1.+redshift)*(1.+redshift)*(1.+redshift)*(1.-simParam->Y);
    }else{
        mean_numdensity_H = simParam->mean_density*(1.+redshift)*(1.+redshift)*(1.+redshift)*(1.-simParam->Y)/(1.-0.75*simParam->Y);
    }
    
    correct_HeII = (1.-0.75*simParam->Y)/(1.-simParam->Y);
    //compute modified photoionization fraction & compute new ionization fraction in ionized regions

    for(int comz=0; comz<local_n0; comz++)
    {
        for(int comy=0; comy<nbins; comy++)
        {
            for(int comx=0; comx<nbins; comx++)
            {
                cell = comz*nbins*nbins + comy*nbins + comx;
                //compute photHI fluctuations (\delta_{photIon})
                photHI = creal(thisGrid->photHI[cell]);
                
                //compute self shielded overdensity
                densSS = ss_calc_densSS(simParam, photHI, temperature, redshift);
                
                //compute modified photHI
                mod_photHI = ss_calc_modPhotHI(creal(thisGrid->igm_density[cell]), densSS);
                
                //compute new XHII
                thisGrid->XHII[cell] = ss_calc_XHII(creal(thisGrid->igm_density[cell])*mean_numdensity_H*correct_HeII*creal(thisGrid->igm_clump[cell]), mod_photHI*creal(thisGrid->photHI[cell]), temperature, simParam->Y) + 0.*I;
            }
        }
    }
}


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* ADDITIONAL FUNCTIONS CURRENTLY NOT USED                                    */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

double calc_mean_photoionization_ionized_field(grid_t *thisGrid)
{
    int nbins = thisGrid->nbins;
    int local_n0 = thisGrid->local_n0;
    
    double sum = 0.;
#ifdef __MPI
    double sum_all = 0.;
#endif
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                sum +=  creal(thisGrid->XHII[i*nbins*nbins+j*nbins+k])*creal(thisGrid->photHI[i*nbins*nbins+j*nbins+k]) + 0.*I;
            }
        }
    }
    
#ifdef __MPI
    MPI_Allreduce(&sum, &sum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sum = sum_all / (double)(nbins*nbins*nbins);
#else
    sum = sum / (double)(nbins*nbins*nbins);
#endif
    
    if(sum<=0.) sum = thisGrid->mean_photHI / (double)(nbins*nbins*nbins);
    
    return sum;
}

double calc_factor_photoionization_ionfraction(grid_t *thisGrid, confObj_t simParam)
{
    int nbins = thisGrid->nbins;
    int local_n0 = thisGrid->local_n0;
    
    double mean_XHII;
    double mean_numdensity_H;
    double chi;
    const double z = simParam->redshift;

    mean_XHII = get_mean_grid(thisGrid->XHII, nbins, local_n0);
        
    //compute mean hydrogen & helium number density
    if(simParam->default_mean_density == 1){
        mean_numdensity_H = 3.*SQR(H0)/(8.*M_PI*G)/mp_g*simParam->omega_b*(1.+z)*(1.+z)*(1.+z)*(1.-simParam->Y);
    }else{
        mean_numdensity_H = simParam->mean_density*(1.+z)*(1.+z)*(1.+z)*(1.-simParam->Y)/(1.-0.75*simParam->Y);
    }
    
    chi = 0.25*simParam->Y/(1.-simParam->Y);

    return mean_numdensity_H*recomb_HII*(1.+chi)*mean_XHII*mean_XHII/(1.-mean_XHII);
}
