#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>    //included because M_PI is not defined in <math.h>
#include <assert.h>
#include <complex.h>
#include <gsl/gsl_integration.h>

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

#include "density_distribution.h"
#include "recombination.h"

#define SQR(X) ((X) * (X))
#define CUB(X) ((X) * (X) * (X))

#ifdef NDEBUG
#define XASSERT(EXP, ...)                                do{} while(0)
#else
#define XASSERT(EXP, ...)                                              \
    do { if (!(EXP)) {                                                  \
            printf("Error in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
            printf(__VA_ARGS__);                                        \
            printf(ANSI_COLOR_BLUE "Bug in code: email Anne Hutter <ahutter@swin.edu.au>"ANSI_COLOR_RESET"\n"); \
            fflush(stdout);                                             \
            exit(EXIT_FAILURE);                                         \
        } \
    } while (0)
#endif
    
    
double get_Xe(double XHII, double XHeII, double XHeIII, double Y, int specie)
{
    double result = 1.;
    
    switch(specie)
    {
        case 0:
            if(XHII > 0.) result = 1. + 0.25*Y/(1.-Y) * (XHeII + 2.*XHeIII)/XHII;
            else result = 1.;
            break;
        case 1:
            if(XHeII > 0.) result = 1. + 4.*(1.-Y)/Y*XHII/XHeII + 2.*XHeIII/XHeII;
            else result = 1.;
            break;
        case 2:
            if(XHeIII > 0.) result = 2. + 4.*(1.-Y)/Y*XHII/XHeIII + XHeII/XHeIII;
            else result = 2.;
            break;
    }
    return result;
}

void compute_number_recombinations(grid_t *thisGrid, confObj_t simParam)
{
    int nbins;
    int local_n0;
    
    double dens;
    double clump;
    double photHI;
    double temp;
    double zstart;
    double redshift;
    
    double Y = simParam->Y;
    double XHII = 1.;
    double XHeII = 1.;
    double XHeIII = 0.;
    double Xe = get_Xe(XHII, XHeII, XHeIII, Y, 0);
    
    nbins = thisGrid->nbins;
    local_n0 = thisGrid->local_n0;
    
    temp = 1.e4;
    if(simParam->calc_ion_history == 1)
    {
        zstart = simParam->redshift_prev_snap;
    }else{
        zstart = 15.;
    }
    redshift = simParam->redshift;
        
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                dens = creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k]);
                clump = creal(thisGrid->igm_clump[i*nbins*nbins+j*nbins+k]);
                photHI = creal(thisGrid->photHI[i*nbins*nbins+j*nbins+k]);

                XHII = creal(thisGrid->XHII[i*nbins*nbins+j*nbins+k]);

                if(simParam->solve_He == 1)
                {
                    XHeIII = creal(thisGrid->XHeIII[i*nbins*nbins+j*nbins+k]);
                    XHeII = creal(thisGrid->XHeII[i*nbins*nbins+j*nbins+k]);
                    
                    Xe = get_Xe(XHII, XHeII, XHeIII, Y, 1);
                    thisGrid->nrec_HeI[i*nbins*nbins+j*nbins+k] =  get_nrec_HeI_history(simParam, dens, clump, 0., temp, zstart, redshift, Xe)+0.*I;
                    
                    Xe = get_Xe(XHII, XHeII, XHeIII, Y, 2);
                    thisGrid->nrec_HeII[i*nbins*nbins+j*nbins+k] =  get_nrec_HeII_history(simParam, dens, clump, 0., temp, zstart, redshift, Xe)+0.*I;
                    
                    Xe = get_Xe(XHII, XHeII, XHeIII, Y, 0);
                }

                if(photHI > 0.)
                {
                    thisGrid->nrec[i*nbins*nbins+j*nbins+k] =  get_nrec_history(simParam, dens, clump, photHI, temp, zstart, redshift, Xe)+0.*I;
                }
                else
                {
                    thisGrid->nrec[i*nbins*nbins+j*nbins+k] = 0.;
                }

#ifdef DEBUG_NREC
                if(simParam->solve_He == 1)
                {
                    if(XHII > 0.) printf("XHII = %e\tXHeII = %e\t XHeIII = %e\tXe = %e\t %e %e %e\n", XHII, XHeII, XHeIII, Xe, creal(thisGrid->nrec[i*nbins*nbins+j*nbins+k]), creal(thisGrid->nrec_HeI[i*nbins*nbins+j*nbins+k]), creal(thisGrid->nrec_HeII[i*nbins*nbins+j*nbins+k]));
                }else{
                   if(XHII > 0.) printf("XHII = %e\tXHeII = %e\t XHeIII = %e\tXe = %e\t %e \n", XHII, XHeII, XHeIII, Xe, creal(thisGrid->nrec[i*nbins*nbins+j*nbins+k]));
                }
                
                printf("nrec = %e\t dens = %e\t photHI = %e\t %e\n", creal(thisGrid->nrec[i*nbins*nbins+j*nbins+k]), dens, photHI, pow(dens, 1./3.));
#endif
                
            }
        }
    }
}

double get_nrec_history(confObj_t simParam, double dens, double clump, double photHI, double temp, double zstart, double redshift, double Xe)
{
    const double time = time_from_redshift_flatuniverse(simParam, redshift, zstart);
    double mean_numdensity_H;
    double factor;
    double factork;
    double XHII;
    double tmp;
    
    if(simParam->default_mean_density == 1){
        mean_numdensity_H = 3.*SQR(H0)/(8.*M_PI*G)/mp_g*simParam->omega_b*CUB(1.+redshift)*(1.-simParam->Y);
    }else{
        mean_numdensity_H = simParam->mean_density*CUB(1.+redshift)*(1.-simParam->Y)/(1.-0.75*simParam->Y);
    }
    
    if(Xe > 0.){
        factor = recomb_HII * Xe * dens * mean_numdensity_H;
        factork = photHI / (clump * factor);
        if(photHI > 0.) XHII = 0.5* factork * (sqrt(4./factork+1.) - 1.);
        else XHII = 1.;
        tmp = XHII * factor * time;
    }else{
        XHII = 0.;
        tmp = 0.;
    }
    
    return tmp / (1. + tmp);
}

double get_nrec_HeI_history(confObj_t simParam, double dens, double clump, double photHeI, double temp, double zstart, double redshift, double Xe)
{
    const double time = time_from_redshift_flatuniverse(simParam, redshift, zstart);
    double mean_numdensity_He;

    if(simParam->default_mean_density == 1){
        mean_numdensity_He = 3.*SQR(H0)/(8.*M_PI*G)/mp_g*simParam->omega_b*CUB(1.+redshift)*0.25*simParam->Y;
    }else{
        mean_numdensity_He = simParam->mean_density*CUB(1.+redshift)*0.25*simParam->Y/(1.-0.75*simParam->Y);
    }
    
    const double factor = recomb_HeII * Xe * dens * mean_numdensity_He;
    const double tmp = factor * time;
    
    return tmp / (1. + tmp);
}

double get_nrec_HeII_history(confObj_t simParam, double dens, double clump, double photHeII, double temp, double zstart, double redshift, double Xe)
{
    const double time = time_from_redshift_flatuniverse(simParam, redshift, zstart);
    double mean_numdensity_He;

    if(simParam->default_mean_density == 1){
        mean_numdensity_He = 3.*SQR(H0)/(8.*M_PI*G)/mp_g*simParam->omega_b*CUB(1.+redshift)*0.25*simParam->Y;
    }else{
        mean_numdensity_He = simParam->mean_density*CUB(1.+redshift)*0.25*simParam->Y/(1.-0.75*simParam->Y);
    }
    
    const double factor = recomb_HeIII * Xe * dens * mean_numdensity_He;
    const double tmp = factor * time;
    
    return tmp / (1. + tmp);
}

void compute_number_recombinations_M2000(grid_t *thisGrid, confObj_t simParam, char *filename, const integral_table_t *thisIntegralTable)
{
    int nbins;
    int local_n0;
    
    double dens;
    double clump;
    double temp;
    double photHI;
    
    double zstart;
    double redshift;
    
    double Y = simParam->Y;
    double XHII = 1.;
    double XHeII = 1.;
    double XHeIII = 0.;
    double Xe = get_Xe(XHII, XHeII, XHeIII, Y, 0);
    
    double *integral_table;
    
    nbins = thisGrid->nbins;
    local_n0 = thisGrid->local_n0;
    
    temp = 1.e4;
    if(simParam->calc_ion_history == 1)
    {
        zstart = simParam->redshift_prev_snap;
    }else{
        zstart = 15.;
    }
    redshift = simParam->redshift;
        
    integral_table = read_table_integrals(filename, thisIntegralTable);
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                dens = creal(thisGrid->igm_density[i*nbins*nbins+j*nbins+k]);
                photHI = creal(thisGrid->photHI[i*nbins*nbins+j*nbins+k]);
                clump = creal(thisGrid->igm_clump[i*nbins*nbins+j*nbins+k]);

                XHII = creal(thisGrid->XHII[i*nbins*nbins+j*nbins+k]);

                if(simParam->solve_He == 1)
                {
                    XHeIII = creal(thisGrid->XHeIII[i*nbins*nbins+j*nbins+k]);
                    XHeII = creal(thisGrid->XHeII[i*nbins*nbins+j*nbins+k]);
                    
                    Xe = get_Xe(XHII, XHeII, XHeIII, Y, 1);
                    thisGrid->nrec_HeI[i*nbins*nbins+j*nbins+k] =  get_nrec_HeI_history(simParam, dens, clump, 0., temp, zstart, redshift, Xe)+0.*I;
                    
                    Xe = get_Xe(XHII, XHeII, XHeIII, Y, 2);
                    thisGrid->nrec_HeII[i*nbins*nbins+j*nbins+k] =  get_nrec_HeII_history(simParam, dens, clump, 0., temp, zstart, redshift, Xe)+0.*I;
                }
                
                if(photHI > 0.)
                {
                    thisGrid->nrec[i*nbins*nbins+j*nbins+k] =  get_nrec_history_M2000(simParam, thisIntegralTable, integral_table, dens, photHI, temp, zstart, redshift)+0.*I;
                }
                else
                {
                    thisGrid->nrec[i*nbins*nbins+j*nbins+k] = 0.;
                }
#ifdef DEBUG_NREC
                printf("nrec = %e\t dens = %e\t photHI = %e\t %e\n", creal(thisGrid->nrec[i*nbins*nbins+j*nbins+k]), dens, photHI, pow(dens, 1./3.));
#endif
                
            }
        }
    }
    
    free(integral_table);
    
//     compute_table_norm_pdf(0., 10., 0.1, 0, 1, "norm_pdf_z0_10.dat");
}

double get_nrec_history_M2000(confObj_t simParam, const integral_table_t *thisIntegralTable, double *integral_table, double dens, double photHI, double temp, double zstart, double redshift)
{
    double tmp;
    double correctFact;
    double factor, dcell;
    int numz, zstart_index, redshift_index;
    int numf, factor_index;
    int numdcell, dcell_index;
    
    if(simParam->default_mean_density == 1){
        correctFact = (1.-0.75*simParam->Y)*(3.*SQR(H0))/(8.*M_PI*G*SQR(simParam->h))/1.8791e-29;

    }else{
        correctFact = (1.-0.75*simParam->Y)*simParam->mean_density/(1.8791e-29/mp_g*simParam->h*simParam->h*simParam->omega_b);
    }

    dcell = dens;
    numdcell = (thisIntegralTable->dcellmax - thisIntegralTable->dcellmin)/thisIntegralTable->ddcell+1;
    dcell_index = (log10(dcell) - thisIntegralTable->dcellmin)/thisIntegralTable->ddcell;
    
    if(dcell_index<0)
    {
        printf("dcell = %e\tdcell_index = %d, not within limits of %d to %d\n", dcell, dcell_index, 0, numdcell);
        dcell_index = 0;
    }
    if(dcell_index>=numdcell)
    {
        printf("dcell = %e\tdcell_index = %d, not within limits of %d to %d\n", dcell, dcell_index, 0, numdcell);
        dcell_index = numdcell-1;
    }
    assert(dcell_index>=0 && dcell_index<numdcell);
    
    factor = (recomb_HII*correctFact)/photHI;
        
    numf = (thisIntegralTable->fmax - thisIntegralTable->fmin)/thisIntegralTable->df+1;
    factor_index = (log10(factor) - thisIntegralTable->fmin)/thisIntegralTable->df;
    
    if(factor_index<0)
    {
        printf("factor = %e\tfactor_index = %d, not within limits of %d to %d\t", factor, factor_index, 0, numf);
        printf(" photHI = %e\n", photHI);
        factor_index = 0;
    }
    if(factor_index>=numf)
    {
        printf("factor = %e\tfactor_index = %d, not within limits of %d to %d\t", factor, factor_index, 0, numf);
        printf(" photHI = %e\n", photHI);
        factor_index = numf-1;
    }
    assert(factor_index>=0 && factor_index<numf);

    numz = (thisIntegralTable->zmax - thisIntegralTable->zmin)/thisIntegralTable->dz+1;
    zstart_index = (zstart  - thisIntegralTable->zmin)/thisIntegralTable->dz;
    redshift_index = (redshift - thisIntegralTable->zmin)/thisIntegralTable->dz;
    
    if(redshift_index<0)
    {
        printf("redshift_index = %d, not within limits of %d to %d\n", redshift_index, 0, numz);
        redshift_index = 0;
    }
    if(redshift_index>=numz)
    {
        printf("redshift_index = %d, not within limits of %d to %d\n", redshift_index, 0, numz);
        redshift_index = numz-1;
    }
    assert(redshift_index>=0 && redshift_index<numz);

    for(int k=0; k<numdcell; k++)
    {
        for(int j=0; j<numf; j++)
        {
            tmp = 0.;
            for(int i=redshift_index; i<zstart_index; i++)
            {
                tmp += integral_table[numz*numf*k + numz*j + i];
            }
            tmp = tmp*photHI*1.e12;
#ifdef DEBUG_NREC
            printf("factor = %d\t dcell_index = %d\t tmp = %e\n",j, k, tmp);
#endif
        }
    }
    
    tmp = 0.;
    for(int i=redshift_index; i<zstart_index; i++)
    {
        tmp += integral_table[numz*numf*dcell_index + numz*factor_index + i];
//         if(photHI>1.e-15) printf("z = %e:\ttmp = %e\n", i*thisIntegralTable->dz + thisIntegralTable->zmin, tmp);
    }
    tmp = tmp*photHI*1.e12;
    
//     printf("dcell = %e\t photHI = %e\t tmp = %e\t zstart = %e\t zend = %e\n", dcell, photHI, tmp, zstart, redshift);
    if(tmp < 0.) tmp = 0.;
    
    return tmp;
}

void compute_number_recombinations_const(grid_t *thisGrid, confObj_t simParam)
{
    int nbins;
    int local_n0;
    
    double zstart;
    double redshift;
    double nrec;
    
    nbins = thisGrid->nbins;
    local_n0 = thisGrid->local_n0;
    
    if(simParam->calc_ion_history == 1)
    {
        zstart = simParam->redshift_prev_snap;
    }else{
        zstart = 15.;
    }
    redshift = simParam->redshift;
    
    // computing recombinations for HII
    nrec = get_nrec_history_constantInTime(simParam, redshift, zstart);
    printf("\n nrec = %e Myr^-1\n", nrec);
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                thisGrid->nrec[i*nbins*nbins+j*nbins+k] =  nrec + 0.*I;
            }
        }
    }

    
    if(simParam->solve_He == 1)
    {
        // computing recombinations for HeII
        nrec = get_nrec_HeI_history_constantInTime(simParam, redshift, zstart);
        printf("\n nrec_HeI = %e Myr^-1\n", nrec);
        for(int i=0; i<local_n0; i++)
        {
            for(int j=0; j<nbins; j++)
            {
                for(int k=0; k<nbins; k++)
                {
                    thisGrid->nrec_HeI[i*nbins*nbins+j*nbins+k] =  nrec + 0.*I;
                }
            }
        }

        // computing recombinations for HeIII
        nrec = get_nrec_HeII_history_constantInTime(simParam, redshift, zstart);
        printf(" nrec_HeII = %e Myr^-1\n", nrec);
        for(int i=0; i<local_n0; i++)
        {
            for(int j=0; j<nbins; j++)
            {
                for(int k=0; k<nbins; k++)
                {
                    thisGrid->nrec_HeII[i*nbins*nbins+j*nbins+k] =  nrec + 0.*I;
                }
            }
        }
    }
}

double get_nrec_history_constantInTime(confObj_t simParam, double z, double zstart)
{
    const double time = time_from_redshift_flatuniverse(simParam, z, zstart)/Myr_s;
    const double nrec = simParam->dnrec_dt * time;
    
    return nrec;
}

double get_nrec_HeI_history_constantInTime(confObj_t simParam, double z, double zstart)
{
    const double time = time_from_redshift_flatuniverse(simParam, z, zstart)/Myr_s;
    const double nrec = simParam->dnrec_HeI_dt * time;
    
    return nrec;
}

double get_nrec_HeII_history_constantInTime(confObj_t simParam, double z, double zstart)
{
    const double time = time_from_redshift_flatuniverse(simParam, z, zstart)/Myr_s;
    const double nrec = simParam->dnrec_HeII_dt * time;
    
    return nrec;
}

//------------------------------------------------------------------------------
// table for pdf
//------------------------------------------------------------------------------

double *read_table_norm_pdf(char *filename)
{
    int len_byte;
    FILE * fp;
    
    int num;
    double *data;
    
    fp = fopen(filename, "rb");
    if(fp == NULL)
    {
        fprintf(stderr, "Error: NO File!");
        exit(EXIT_FAILURE);
    }
    fseek(fp, 0, SEEK_END);
    len_byte = ftell(fp);
    rewind(fp);
    
    num = len_byte/(sizeof(double)*3);
    data = malloc(num*3*sizeof(double));
    if(data == NULL)
    {
        fprintf(stderr, "data in read_table_norm_pdf (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    if(fread(data, sizeof(double), 3*num, fp) != (unsigned int)num)
    {
        fprintf(stderr, "Error: Could not read pdf table!\n");
        exit(EXIT_FAILURE);
    }
    fclose(fp);
    
    fp = fopen("norm_pdf_txt.dat", "w");
    for(int i=0; i<num; i++)
    {
        fprintf(fp, "%e\t%e\t%e\n", data[3*i], data[3*i+1], data[3*i+2]);
    }
    fclose(fp);
    
    return data;
}

void compute_table_norm_pdf(double zmin, double zmax, double d, int rank, int size, char *filename)
{
    pdf_params_t *pdf_params;
    double *array;
    int num, num_size, offset;
    int offset_write;
    double redshift;
    
    num = (zmax - zmin)/d;
    int tmp = num/size;
    
    if(size*tmp == num)
    {
        num_size = tmp;
        offset = 0;
    }else
    {
        num_size = ((num%size)>rank) ? tmp+1 : tmp;
        offset = ((num%size)>rank) ? rank*num_size -rank+num%size : rank*num_size;
    }
    
    offset_write = num_size*rank + offset;
    
    pdf_params = malloc(sizeof(pdf_params_t));
    if(pdf_params == NULL)
    {
        fprintf(stderr, "pdf_params in compute_table_norm_pdf (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    array = malloc(4*num_size*sizeof(double));
    if(array == NULL)
    {
        fprintf(stderr, "array in compute_table_norm_pdf (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i=0; i<num_size; i++)
    {
        redshift = zmin + d*(i + offset_write);
        dd_set_norm_pdf(pdf_params, redshift);
        array[4*i] = redshift;
        array[4*i+1] = pdf_params->amplitude;
        array[4*i+2] = pdf_params->constant;
        array[4*i+3] = pdf_params->beta;
    }
    
#ifdef __MPI
    write_table(num_size*4, offset_write*4, array, filename);
#else
    write_table(num_size*4, array, filename);
#endif
    
    FILE *f;
    f=fopen(filename, "w");
    for(int i=0; i<num_size; i++)
    {
        fprintf(f,"%e\t%e\t%e\t%e\n",array[4*i],array[4*i+1],array[4*i+2],array[4*i+3]);
    }
    fclose(f);
    
    free(array);
}

double amplitude_norm_pdf(double z)
{
    return (0.053*z+0.03)*(1.+exp(-0.72*z+1.9));
}

double constant_norm_pdf(double z)
{
    return (1.-exp(-0.66*z+1.9));
}

//------------------------------------------------------------------------------
// table for integrals over density and redshift
//------------------------------------------------------------------------------

double *read_table_integrals(char *filename, const integral_table_t *thisIntegralTable)
{
    int len_byte;
    FILE * fp;
    
    int num;
    double *array;
    int numz, numf, numdcell;

    fp = fopen(filename, "rb");
    if(fp == NULL)
    {
        printf("Error: NO File!");
        exit(1);
    }
    fseek(fp, 0, SEEK_END);
    len_byte = ftell(fp);
    rewind(fp);
    
    num = len_byte/(sizeof(double));
    array = malloc(num*sizeof(double));
    if(array == NULL)
    {
        fprintf(stderr, "array in read_table_integrals (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    if(fread(array, sizeof(double), num, fp) != (unsigned int)num)
    {
        fprintf(stderr, "Error: Could not read nrec integral table!\n");
        exit(EXIT_FAILURE);
    }
    fclose(fp);
    
    numz = (thisIntegralTable->zmax - thisIntegralTable->zmin)/thisIntegralTable->dz+1;
    
    numf = (thisIntegralTable->fmax - thisIntegralTable->fmin)/thisIntegralTable->df+1;
    
    numdcell = (thisIntegralTable->dcellmax - thisIntegralTable->dcellmin)/thisIntegralTable->ddcell+1;
    
    assert(num == numz*numf*numdcell);
    
//     for(int i=0; i<numdcell; i++) printf("%e\n",array[i]);
    
    return array;
}

integral_table_t * initIntegralTable(double zmin, double zmax, double dz, double fmin, double fmax, double df, double dcellmin, double dcellmax, double ddcell)
{
    integral_table_t *newIntegralTable;
    newIntegralTable = malloc(sizeof(integral_table_t));
    if(newIntegralTable == NULL)
    {
        fprintf(stderr, "newIntegralTable in initIntegralTable (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    newIntegralTable->zmin = zmin;
    newIntegralTable->zmax = zmax;
    newIntegralTable->dz = dz;
    
    newIntegralTable->fmin = fmin;
    newIntegralTable->fmax = fmax;
    newIntegralTable->df = df;
    
    newIntegralTable->dcellmin = dcellmin;
    newIntegralTable->dcellmax = dcellmax;
    newIntegralTable->ddcell = ddcell;
    
    return newIntegralTable;
}

//------------------------------------------------------------------------------
// table to perform integral over density
//------------------------------------------------------------------------------

dens_table_t *initDensTable(double constant_min, double constant_max, double dconstant, double factor_min, double factor_max, double dfactor)
{
    dens_table_t *newDensTable;
    newDensTable = malloc(sizeof(dens_table_t));
    if(newDensTable == NULL)
    {
        fprintf(stderr, "newDensTable in initDensTable (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    newDensTable->constant_min = constant_min;
    newDensTable->constant_max = constant_max;
    newDensTable->dconstant = dconstant;
    
    newDensTable->factor_min = factor_min;
    newDensTable->factor_max = factor_max;
    newDensTable->dfactor = dfactor;
    
    return newDensTable;
}

double dens_integrand(double x, void * p)
{
    dens_integrand_t * params = (dens_integrand_t *)p;
    
    double frac = 2./3;
    double constant = params->constant;
    double factor = params->factor;
    
    double nrec = (sqrt(1.+4.*x*factor)-1.);
    
    return exp(-SQR(pow(x,-frac)-constant))*pow(x,-2.5)*nrec;
}

double calc_dens_integral(dens_integrand_t *params)
{
    gsl_function F;
    F.function = &dens_integrand;
    F.params = (void *)params;
    double result, error;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
    
    gsl_integration_qagiu(&F, 0., 1.e-9, 1.e-9, 10000, w, &result, &error);
    
    gsl_integration_workspace_free(w);

    return result;
}

double *create_table_dens(dens_table_t *thisDensTable)
{
    double *array;
    dens_integrand_t *params;
    
    double constant, factor;
    double dconstant, dfactor;
    int num1, num2;
    
    dconstant = thisDensTable->dconstant;
    num1 = (thisDensTable->constant_max - thisDensTable->constant_min)/dconstant;
    
    dfactor = thisDensTable->dfactor;
    num2 = (thisDensTable->factor_max - thisDensTable->factor_min)/dfactor;
    
    array = malloc(num1*num2*3*sizeof(double));
    if(array == NULL)
    {
        fprintf(stderr, "array in create_table_dens (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    params = malloc(sizeof(dens_integrand_t));
    if(params == NULL)
    {
        fprintf(stderr, "params in create_table_dens (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i=0; i<num1; i++)
    {
        constant = thisDensTable->constant_min+dconstant*i;
        params->constant = constant;
        for(int j=0; j<num2; j++)
        {
            factor = pow(10.,thisDensTable->factor_min+dfactor*j);
            params->factor = factor;
            array[i*num2*3+3*j] = constant;
            array[i*num2*3+3*j+1] = factor;
            array[i*num2*3+3*j+2]= calc_dens_integral(params);
        }
    }
    
    free(params);
    
    return array;
}

void compute_table_dens(double constant_min, double constant_max, double d1, double factor_min, double factor_max, double d2, int rank, int size, char *filename)
{
    double *array;
    dens_integrand_t *params;
    
    double constant, factor;
    int num1, num1_size, num2, offset, tmp;
    
    int offset_write;
    
    d1 = 0.001;
    num1 = (constant_max - constant_min)/d1;
    
    d2 = 0.01;
    num2 = (factor_max - factor_min)/d2;
    
    tmp = num1/size;
    
    if(size*tmp == num1)
    {
        num1_size = tmp;
        offset = 0;
    }else
    {
        num1_size = ((num1%size)>rank) ? tmp+1 : tmp;
        offset = ((num1%size)>rank) ? rank*num1_size -rank+num1%size : rank*num1_size;
    }
    
    offset_write = num1_size*rank + offset;
    
    array = malloc(num1*num2*3*sizeof(double));
    if(array == NULL)
    {
        fprintf(stderr, "array in compute_table_dens (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    params = malloc(sizeof(dens_integrand_t));
    if(params == NULL)
    {
        fprintf(stderr, "params in compute_table_dens (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i=0; i<num1_size; i++)
    {
        constant = constant_min+d1*(i + offset_write);
        params->constant = constant;
        for(int j=0; j<num2; j++)
        {
            factor = pow(10.,factor_min+d2*j);
            params->factor = factor;
            array[i*num2*3+3*j] = constant;
            array[i*num2*3+3*j+1] = factor;
            array[i*num2*3+3*j+2]= calc_dens_integral(params);
        }
    }

#ifdef __MPI
    write_table(num1*num2*3, offset_write*num2*3, array, filename);
#else
    write_table(num1*num2*3, array, filename);
#endif
    
    free(params);
    free(array);
}

//------------------------------------------------------------------------------
// table to perform integral over redshift
//------------------------------------------------------------------------------

redshift_table_t *initRedshiftTable(double dens_cell_min, double dens_cell_max, double ddens_cell, double zmin, double zmax, double dz)
{
    redshift_table_t *newRedshiftTable;
    newRedshiftTable = malloc(sizeof(redshift_table_t));
    if(newRedshiftTable == NULL)
    {
        fprintf(stderr, "newRedshiftTable in initRedshiftTable (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    newRedshiftTable->dens_cell_min = dens_cell_min;
    newRedshiftTable->dens_cell_max = dens_cell_max;
    newRedshiftTable->ddens_cell = ddens_cell;
    
    newRedshiftTable->zmin = zmin;
    newRedshiftTable->zmax = zmax;
    newRedshiftTable->dz = dz;
    
    return newRedshiftTable;
}

double redshift_integrand(double x, void * p)
{
    redshift_integrand_t * params = (redshift_integrand_t *)p;
    
    double frac = 2./3.;
    double dens_cell = params->dens_cell;
    confObj_t simParam = params->simParam;

    double tmp = H0*sqrt((simParam->omega_b*CUB(1.+x)+simParam->omega_l))*(1.+x);
    
//     printf("tmp = %e\t%e\t%e\t%e\n", tmp, dens_cell, x, 2.*SQR(frac*7.61/((1.+x)*dens_cell)));
    return exp(2.*SQR(frac*7.61/((1.+x)*dens_cell)))*amplitude_norm_pdf(x)/tmp;
}

double calc_redshift_integral(redshift_integrand_t *params, double d2)
{
    gsl_function F;
    F.function = &redshift_integrand;
    F.params = (void *)params;
    double result, error;
    double z = params->z;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000); 
    
    gsl_integration_qag(&F, z, z+d2, 1.e-5, 1.e-3, 10000, 1, w, &result, &error);
    
    gsl_integration_workspace_free(w);

    return result;
}

double *create_table_redshift(redshift_table_t *thisRedshiftTable, confObj_t simParam)
{
    double * array;
    redshift_integrand_t *params;
    
    double dens_cell, z;
    double ddens_cell, dz;
    int num1, num2;
    
    ddens_cell = thisRedshiftTable->ddens_cell;
    num1 = (thisRedshiftTable->dens_cell_max - thisRedshiftTable->dens_cell_min)/ddens_cell;

    dz = thisRedshiftTable->dz;
    num2 = (thisRedshiftTable->zmax - thisRedshiftTable->zmin)/dz;
    
    array = malloc(num1*num2*3*sizeof(double));
    if(array == NULL)
    {
        fprintf(stderr, "array in create_table_redshift (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    params = malloc(sizeof(redshift_integrand_t));
    if(params == NULL)
    {
        fprintf(stderr, "params in create_table_redshift (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    params->simParam = simParam;
    
    for(int i=0; i<num1; i++)
    {
        dens_cell = pow(10.,thisRedshiftTable->dens_cell_min+thisRedshiftTable->ddens_cell*i);
        params->dens_cell = dens_cell;
        for(int j=0; j<num2; j++)
        {
            z = thisRedshiftTable->zmin+j*dz;
            params->z = z;
            array[i*num2*3+3*j] = dens_cell;
            array[i*num2*3+3*j+1] = z;
            array[i*num2*3+3*j+2] = calc_redshift_integral(params, dz);
        }
    }
    
    free(params);
    
    return array;
}

void compute_table_redshift(double dens_cell_min, double dens_cell_max, double d1, double zmin, double zmax, double d2, confObj_t simParam, int rank, int size, char *filename)
{
    double * array;
    redshift_integrand_t *params;
    
    double dens_cell, z;
    int num1, num1_size, num2, offset, tmp;
    
    int offset_write;
    
    d1 = 0.01;
    num1 = (dens_cell_max - dens_cell_min)/d1;

    d2 = 0.01;
    num2 = (zmax - zmin)/d2;
    
    tmp = num1/size;
    
    if(size*tmp == num1)
    {
        num1_size = tmp;
        offset = 0;
    }else
    {
        num1_size = ((num1%size)>rank) ? tmp+1 : tmp;
        offset = ((num1%size)>rank) ? rank*num1_size -rank+num1%size : rank*num1_size;
    }
    
    offset_write = num1_size * rank + offset;
    
    array = malloc(num1*num2*3*sizeof(double));
    if(array == NULL)
    {
        fprintf(stderr, "array in compute_table_redshift (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    params = malloc(sizeof(redshift_integrand_t));
    if(params == NULL)
    {
        fprintf(stderr, "params in compute_table_redshift (recombination.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    params->simParam = simParam;
    
    for(int i=0; i<num1_size; i++)
    {
        dens_cell = pow(10.,dens_cell_min+d1*(i + offset_write));
        params->dens_cell = dens_cell;
        for(int j=0; j<num2; j++)
        {
            z = zmin+j*d2;
            params->z = z;
            array[i*num2*3+3*j] = dens_cell;
            array[i*num2*3+3*j+1] = z;
            array[i*num2*3+3*j+2] = calc_redshift_integral(params, d2);
        }
    }

#ifdef __MPI
    write_table(num1*num2*3, offset_write*num2*3, array, filename);
#else
    write_table(num1*num2*3, array, filename);
#endif
    
    free(params);
    free(array);
}

//------------------------------------------------------------------------------
// write table
//------------------------------------------------------------------------------

#ifdef __MPI
void write_table(int num, int offset, double *array, char *filename)
{
    MPI_File mpifile;
    MPI_Offset offset_mpi;
    MPI_Status status;
    
    offset_mpi = offset*sizeof(double);
    
    printf("num = %d\t offset = %d\n", num, offset);
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifile);
    MPI_File_write_at_all(mpifile, offset_mpi, array, num, MPI_DOUBLE, &status);
    MPI_File_close(&mpifile);
}
#else
void write_table(int num, double *array, char *filename)
{
    FILE * fp;
    
    fp = fopen(filename, "wb");
    fwrite(array, sizeof(double), num, fp);
    fclose(fp);
}
#endif
