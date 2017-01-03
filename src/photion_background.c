#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>	//included because M_PI is not defined in <math.h>
#include <assert.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "photion_background.h"
#include "phys_const.h"
#include "confObj.h"
#include "grid.h"
#include "sources.h"


photIonlist_t *allocate_photIonlist(int Nallocated)
{
    photIonlist_t *newPhotIonlist;
	newPhotIonlist = malloc(sizeof(photIonlist_t));
	if(newPhotIonlist == NULL)
	{
		fprintf(stderr, "newPhotIonlist in allocate_photIonlist (photion_background.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	newPhotIonlist->num = 0;
	newPhotIonlist->Nallocated = Nallocated;
	newPhotIonlist->redshift = NULL;
    newPhotIonlist->photHI = NULL;
    newPhotIonlist->QHII = NULL;
    
    return newPhotIonlist;
}

void deallocate_photIonlist(photIonlist_t *thisPhotIonlist)
{
    if(thisPhotIonlist != NULL)
    {
        if(thisPhotIonlist->redshift != NULL) free(thisPhotIonlist->redshift);
        if(thisPhotIonlist->photHI != NULL) free(thisPhotIonlist->photHI);
        if(thisPhotIonlist->QHII != NULL) free(thisPhotIonlist->QHII);
        free(thisPhotIonlist);
    }
}

photIonlist_t *read_photIonlist(char *filename)
{
    FILE * fp;
    char line[256];
    int i;
    
    int numLines;
    photIonlist_t *newPhotIonlist;
    
    double redshift, photionHI, photheatHI, QHII;
        
    fp = fopen(filename, "rt");
    if(fp == NULL)
    {
        fprintf(stderr, " Can not open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
    
    if(fgets(line, 256, fp) != NULL)
    {
        sscanf(line,"%d\n", &numLines);
    }
    
    newPhotIonlist = allocate_photIonlist(numLines);
    newPhotIonlist->num = numLines;
    newPhotIonlist->Nallocated = numLines;
    
    newPhotIonlist->redshift = malloc(numLines * sizeof(double));
    newPhotIonlist->photHI = malloc(numLines * sizeof(double));
    newPhotIonlist->QHII = malloc(numLines * sizeof(double));

    i = 0;
    while(fgets(line,256,fp)!=NULL)
	{
        sscanf(line,"%lf\t%lf\t%lf\t%lf\n", &redshift, &photionHI, &photheatHI, &QHII);
        newPhotIonlist->redshift[i] = redshift;
        newPhotIonlist->photHI[i] = photionHI;
        newPhotIonlist->QHII[i] = QHII;
        i++;
    }
    
    fclose(fp);
    
    return newPhotIonlist;
}

double get_photHI_from_redshift(photIonlist_t *thisPhotIonlist, double redshift)
{
    int numLines = thisPhotIonlist->num;
    int i;
    double value = 0.;
    double zmin, zmax, photIonMin, photIonMax;
    
    for(i=0; i< numLines; i++)
    {
        if(thisPhotIonlist->redshift[i] >= redshift) break;
    }
    if(i < numLines-1)
    {
        zmin = thisPhotIonlist->redshift[i];
        zmax = thisPhotIonlist->redshift[i+1];
        photIonMin = thisPhotIonlist->photHI[i];
        photIonMax = thisPhotIonlist->photHI[i+1];
        value = photIonMin + (photIonMax - photIonMin)/(zmax - zmin) * (redshift - zmin);
    }else{
        value = thisPhotIonlist->photHI[i];
    }
    
    return value;
}

double get_photHI_from_fillingfactor(photIonlist_t *thisPhotIonlist, double QHII)
{
    int numLines = thisPhotIonlist->num;
    int i;
    double value = 0.;
    double zmin, zmax, Qmin, Qmax;

    for(i=0; i< numLines; i++)
    {
        if(thisPhotIonlist->QHII[i] <= QHII) break;
    }
    
    if(i > 0)
    {
        zmin = thisPhotIonlist->redshift[i-1];
        zmax = thisPhotIonlist->redshift[i];
        Qmin = thisPhotIonlist->QHII[i-1];
        Qmax = thisPhotIonlist->QHII[i];
        value = Qmin + (Qmax - Qmin)/(zmax - zmin) * (QHII - zmin);
    }else{
        value = thisPhotIonlist->QHII[i];
    }
    
    return value;
}
