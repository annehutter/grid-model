#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>    //included because M_PI is not defined in <math.h>
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

#include "redshift_tools.h"

#define SQR(X) ((X) * (X))

void compute_cum_values(grid_t *thisGrid, confObj_t simParam, int specie)
{
    int nbins;
    int local_n0;
    double box_size;
    
    double evol_time;
    double evol_time_fromPrevSnap;
    double z; 
    double mean_numdensity_H;
    double mean_numdensity_He;
    double h;
    
    double Nion;
    
    nbins = thisGrid->nbins;
    local_n0 = thisGrid->local_n0;
    box_size = thisGrid->box_size;
    
    //compute time from previous snapshot, i.e. the time span to evolve ionization fields
    if(simParam->calc_ion_history == 1)
    {
        evol_time_fromPrevSnap = time_from_redshift_flatuniverse(simParam, simParam->redshift, simParam->redshift_prev_snap);
        
        if(specie != 1 && specie != 2)
        {
            printf("\n zstart = %e\t zend = %e\t evol_time = %e + %e Myrs", simParam->redshift_prev_snap, simParam->redshift, simParam->evol_time, time_from_redshift_flatuniverse(simParam, simParam->redshift, simParam->redshift_prev_snap)/Myr_s);
            evol_time = simParam->evol_time*Myr_s + evol_time_fromPrevSnap;
            simParam->evol_time = evol_time/Myr_s;
        }else{
            evol_time = simParam->evol_time*Myr_s;
            printf("\n evol_time = %e Myrs", evol_time/Myr_s);
        }
    }else{
        evol_time = simParam->evol_time*Myr_s;
        evol_time_fromPrevSnap = evol_time;
        printf("\n evol_time = %e Myrs", evol_time/Myr_s);
    }
    
    z = simParam->redshift;
    
    //compute mean hydrogen & helium number density
    if(simParam->default_mean_density == 1){
        mean_numdensity_H = 3.*SQR(H0)/(8.*M_PI*G)/mp_g*simParam->omega_b*(1.+z)*(1.+z)*(1.+z)*(1.-simParam->Y);
        mean_numdensity_He = 3.*SQR(H0)/(8.*M_PI*G)/mp_g*simParam->omega_b*(1.+z)*(1.+z)*(1.+z)*0.25*simParam->Y;
    }else{
        mean_numdensity_H = simParam->mean_density*(1.+z)*(1.+z)*(1.+z)*(1.-simParam->Y)/(1.-0.75*simParam->Y);
        mean_numdensity_He = simParam->mean_density*(1.+z)*(1.+z)*(1.+z)*0.25*simParam->Y/(1.-0.75*simParam->Y);
    }

    //compute number of ionizing photons and absorptions in each cell
    h = simParam->h;
    const double volume = pow(box_size/(h*(double)nbins*(1.+z))*Mpc_cm,3);
    
    if(specie == 1){
        printf("\n mean_numdensity_He at z=%e is %e cm^-3\n evol_time_fromPrevSnap = %e Myr\n", z, mean_numdensity_He, evol_time_fromPrevSnap/Myr_s);

        for(int i=0; i<local_n0; i++)
        {
            for(int j=0; j<nbins; j++)
            {
                for(int k=0; k<nbins; k++)
                {
                    Nion = creal(thisGrid->nion_HeI[i*nbins*nbins+j*nbins+k])*evol_time_fromPrevSnap;

                    thisGrid->cum_nion_HeI[i*nbins*nbins+j*nbins+k] += Nion + 0.*I;

                    thisGrid->cum_nrec_HeI[i*nbins*nbins+j*nbins+k] += creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity_He*volume*creal(thisGrid->nrec_HeI[i*nbins*nbins+j*nbins+k]);
                    
                    thisGrid->cum_nabs_HeI[i*nbins*nbins+j*nbins+k] = creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity_He*volume + creal(thisGrid->cum_nrec_HeI[i*nbins*nbins+j*nbins+k]);
                }
            }
        }
    }else if(specie == 2){
        printf("\n mean_numdensity_He at z=%e is %e cm^-3\n evol_time_fromPrevSnap = %e Myr\n", z, mean_numdensity_He, evol_time_fromPrevSnap/Myr_s);

        for(int i=0; i<local_n0; i++)
        {
            for(int j=0; j<nbins; j++)
            {
                for(int k=0; k<nbins; k++)
                {
                    Nion = creal(thisGrid->nion_HeII[i*nbins*nbins+j*nbins+k])*evol_time_fromPrevSnap;

                    thisGrid->cum_nion_HeII[i*nbins*nbins+j*nbins+k] += Nion + 0.*I;

                    thisGrid->cum_nrec_HeII[i*nbins*nbins+j*nbins+k] += creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity_He*volume*creal(thisGrid->nrec_HeII[i*nbins*nbins+j*nbins+k]);
                    
                    thisGrid->cum_nabs_HeII[i*nbins*nbins+j*nbins+k] = creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity_He*volume + creal(thisGrid->cum_nrec_HeII[i*nbins*nbins+j*nbins+k]);
                }
            }
        }
    }else{
        printf("\n mean_numdensity_H at z=%e is %e cm^-3\n evol_time_fromPrevSnap = %e Myr\n", z, mean_numdensity_H, evol_time_fromPrevSnap/Myr_s);

        for(int i=0; i<local_n0; i++)
        {
            for(int j=0; j<nbins; j++)
            {
                for(int k=0; k<nbins; k++)
                {
                    Nion = creal(thisGrid->nion[i*nbins*nbins+j*nbins+k])*evol_time_fromPrevSnap;
//                     if(Nion>0) printf("%d %d %d: Nion = %e\n", i, j, k, Nion);

                    thisGrid->cum_nion[i*nbins*nbins+j*nbins+k] += Nion + 0.*I;

                    thisGrid->cum_nrec[i*nbins*nbins+j*nbins+k] += creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity_H*volume*creal(thisGrid->nrec[i*nbins*nbins+j*nbins+k]);
                    
                    thisGrid->cum_nabs[i*nbins*nbins+j*nbins+k] = creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k])*mean_numdensity_H*volume + creal(thisGrid->cum_nrec[i*nbins*nbins+j*nbins+k]);

//                     if(creal(thisGrid->cum_nion[i*nbins*nbins+j*nbins+k])>0. && k<1) printf("Nion = %e\t cumNion = %e\t cumNabs = %e\t volume = %e\t dens = %e\n", Nion, creal(thisGrid->cum_nion[i*nbins*nbins+j*nbins+k]), creal(thisGrid->cum_nabs[i*nbins*nbins+j*nbins+k]), volume, creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k]));

                }
            }
        }
    }
}

// versatile function (thisGrid is only used for grid dimensions and domain decomposition
void compute_Q(grid_t *thisGrid, fftw_complex *frac_Q, fftw_complex *nion, fftw_complex *nabs)
{
    int nbins;
    int local_n0;
    
    nbins = thisGrid->nbins;
    local_n0 = thisGrid->local_n0;
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                frac_Q[i*nbins*nbins+j*nbins+k] = creal(nion[i*nbins*nbins+j*nbins+k])/creal(nabs[i*nbins*nbins+j*nbins+k]) + 0.*I;
            }
        }
    }
}
