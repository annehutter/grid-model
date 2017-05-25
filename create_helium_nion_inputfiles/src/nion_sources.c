#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>    //included because M_PI is not defined in <math.h>
#include <assert.h>

#include "phys_const.h"
#include "utils.h"
#include "sources.h" 
#include "confObj.h"
#include "cross_sections.h"
#include "nion_sources.h"
#include "recombination_coeff.h"

void modify_nion_fesc(sourcelist_t *thisSourcelist, double factorNion, double factorFesc)
{
    int numSources = thisSourcelist->numSources;
    source_t *thisSource;

    for(int i=0; i<numSources; i++)
    {
        thisSource = &(thisSourcelist->source[i]);
                
        thisSource->Nion = factorNion*thisSource->Nion;
        thisSource->fesc = factorFesc*thisSource->fesc;
    }

}

void compute_HeII_nion(confObj_t simParam, sourcelist_t *thisSourcelist, sourcelist_t *thisSourcelist_HeI)
{
    int numSources = thisSourcelist->numSources;
    source_t *thisSource;
    source_t *thisSource_HeI;
    
    double Y = simParam->Y;
    double T = simParam->temperature;
    double spectral_index = simParam->source_slope_index;
    double y = 1./(1. + 0.25*Y/(1.-Y)*crossSecHeI(nu_HeI)/crossSecHI(nu_HeI));
    
//     double factor_HeII = (1.-y) * recomb_HII_B / (recomb_HeII_B + y * recomb_HeII_1) * pow(nu_HeI/nu_HI, -spectral_index);
    double factor_HeII = (1.-y) * recB_HII(T) / (recB_HeII(T) + y * (rec_HeII(T) - recB_HeII(T))) * pow(nu_HeI/nu_HI, -spectral_index);
    
    for(int i=0; i<numSources; i++)
    {
        thisSource = &(thisSourcelist->source[i]);
        
        thisSource_HeI = &(thisSourcelist_HeI->source[i]);
        
        thisSource_HeI->Nion = factor_HeII*thisSource->Nion;
        thisSource_HeI->fesc = thisSource->fesc;
        thisSource_HeI->pos[0] = thisSource->pos[0];
        thisSource_HeI->pos[1] = thisSource->pos[1];
        thisSource_HeI->pos[2] = thisSource->pos[2];
    }

}


void compute_HeIII_nion(confObj_t simParam, sourcelist_t *thisSourcelist, sourcelist_t *thisSourcelist_HeII)
{
    int numSources = thisSourcelist->numSources;
    source_t *thisSource;
    source_t *thisSource_HeII;
    
    double T = simParam->temperature;
    double spectral_index = simParam->source_slope_index;
    
//     double factor_HeIII= recomb_HII_B / recomb_HeIII_B * pow(nu_HeII/nu_HI, -spectral_index);
    double factor_HeIII = recB_HII(T) / recB_HeIII(T) * pow(nu_HeII/nu_HI, -spectral_index);    
    
    for(int i=0; i<numSources; i++)
    {
        thisSource = &(thisSourcelist->source[i]);
        
        thisSource_HeII = &(thisSourcelist_HeII->source[i]);
        
        thisSource_HeII->Nion = factor_HeIII*thisSource->Nion;
        thisSource_HeII->fesc = thisSource->fesc;
        thisSource_HeII->pos[0] = thisSource->pos[0];
        thisSource_HeII->pos[1] = thisSource->pos[1];
        thisSource_HeII->pos[2] = thisSource->pos[2];
    }

}
