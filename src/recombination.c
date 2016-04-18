#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>	//included because M_PI is not defined in <math.h>
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

#include "density_distribution.h"
#include "recombination.h"

#define SQR(X) ((X) * (X))
#define CUB(X) ((X) * (X) * (X))


void compute_number_recombinations(grid_t *thisGrid, confObj_t simParam, char *filename, integral_table_t *thisIntegralTable)
{
  	int nbins;
	int local_n0;
	
	double dens;
	double photHI;
	double temp;
	double zstart;
	double redshift;
	
	double *integral_table;
	
	nbins = thisGrid->nbins;
	local_n0 = thisGrid->local_n0;
	
	temp = 1.e4;
	zstart = 15.;
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
				thisGrid->nrec[i*nbins*nbins+j*nbins+k] =  get_nrec_history(simParam, thisIntegralTable, integral_table, dens, photHI, temp, zstart, redshift)+0.*I;
// 				printf("nrec = %e\t dens = %e\t photHI = %e\t %e\n", creal(thisGrid->nrec[i*nbins*nbins+j*nbins+k]), dens, photHI, pow(dens, 1./3.));
			}
		}
	}
	
	free(integral_table);
}

// double get_nrec_history(confObj_t simParam, double *norm_pdf, dens_table_t *thisDensTable, double *dens_table, redshift_table_t *thisRedshiftTable, double *redshift_table, double dens, double photHI, double temp, double zstart, double redshift)
// {
// 	double tmp, tmp_dens, tmp_redshift;
// 	double mean_numdensity, factor;
// 	double z;
// 	int num_z, zstart_index, redshift_index;
// 	int constant_index;
// 	int num_factor, factor_index;
// 	int dens_cell_index;
// 	
// 	dens_integrand_t *dens_params;
// 	redshift_integrand_t *redshift_params;
// 	
// 	dens_params = malloc(sizeof(dens_integrand_t));
// 	redshift_params = malloc(sizeof(redshift_integrand_t));
// 	
// 	if(simParam->default_mean_density == 1){
// 		mean_numdensity = rho_g_cm/mp_g*simParam->h*simParam->h*simParam->omega_b;
// 	}else{
// 		mean_numdensity = simParam->mean_density;
// 	}
// 	factor = mean_numdensity*recomb_HII/photHI;
// 	printf("numdens = %e\tfactor = %e\n",mean_numdensity, factor);
// 
// 	redshift_params->dens_cell = pow(dens, 1./3.);
// 	redshift_params->simParam = simParam;
// 	
// 	dens_cell_index = (redshift_params->dens_cell - thisRedshiftTable->dens_cell_min)/thisRedshiftTable->ddens_cell;
// 	
// 	num_factor = (thisDensTable->factor_max - thisDensTable->factor_min)/thisDensTable->dfactor;
// 	
// 	num_z = (thisRedshiftTable->zmax - thisRedshiftTable->zmin)/thisRedshiftTable->dz;
// 	zstart_index = (zstart  - thisRedshiftTable->zmin)/thisRedshiftTable->dz;
// 	redshift_index = (redshift - thisRedshiftTable->zmin)/thisRedshiftTable->dz;
// 	
// 	tmp = 0.;
// 	for(int i=redshift_index; i<zstart_index; i++)
// 	{
// 		z = thisRedshiftTable->zmin + i*thisRedshiftTable->dz;
// 		
// 		dens_params->factor = factor*CUB(1.+z);
// 		dens_params->constant = constant_norm_pdf(z);
// 		
// 		redshift_params->z = z;
// 		
// 		factor_index = (log10(dens_params->factor) - thisDensTable->factor_min)/thisDensTable->dfactor;
// 		constant_index = (dens_params->constant - thisDensTable->constant_min)/thisDensTable->dconstant;
// 		
// 		tmp_dens = dens_table[3*num_factor*constant_index+3*factor_index+2];
// 		tmp_redshift = redshift_table[3*num_z*dens_cell_index+3*i+2];
// 		printf("%d: z = %e:\ttmp_dens = %e\t tmp_redshift = %e\n",i, z, tmp_dens, tmp_redshift);
// 		tmp += tmp_dens*tmp_redshift;
// 	}
// 	tmp = tmp*0.5*photHI;
// 	
// 	free(dens_params);
// 	free(redshift_params);
// 	
// 	return tmp;
// }

double get_nrec_history(confObj_t simParam, integral_table_t *thisIntegralTable, double *integral_table, double dens, double photHI, double temp, double zstart, double redshift)
{
	double tmp;
	double mean_numdensity, correctFact;
	double z, factor, dcell;
	int numz, zstart_index, redshift_index;
	int numf, factor_index;
	int numdcell, dcell_index;
	
	if(simParam->default_mean_density == 1){
		mean_numdensity = rho_g_cm/mp_g*simParam->h*simParam->h*simParam->omega_b;
		correctFact = 1.;
	}else{
		mean_numdensity = simParam->mean_density;
		correctFact = simParam->mean_density/(rho_g_cm/mp_g*simParam->h*simParam->h*simParam->omega_b);
	}

	dcell = dens;
	numdcell = (thisIntegralTable->dcellmax - thisIntegralTable->dcellmin)/thisIntegralTable->ddcell+1;
	dcell_index = (log10(dcell) - thisIntegralTable->dcellmin)/thisIntegralTable->ddcell;
	
	assert(dcell_index>=0 && dcell_index<numdcell);
	
	factor = (recomb_HII*correctFact)/photHI;
	numf = (thisIntegralTable->fmax - thisIntegralTable->fmin)/thisIntegralTable->df+1;
	factor_index = (log10(factor) - thisIntegralTable->fmin)/thisIntegralTable->df;
	
	assert(factor_index>=0 && factor_index<numf);

	numz = (thisIntegralTable->zmax - thisIntegralTable->zmin)/thisIntegralTable->dz+1;
	zstart_index = (zstart  - thisIntegralTable->zmin)/thisIntegralTable->dz;
	redshift_index = (redshift - thisIntegralTable->zmin)/thisIntegralTable->dz;
	
	assert(redshift_index>=0 && redshift_index<numz);

// 	printf("numdens = %e\tfactor = %e\n",mean_numdensity, factor);
// 	printf("dcell_index = %d\t factor_index = %d\t redshift_index = %d\n", dcell_index, factor_index, redshift_index);

	tmp = 0.;
	for(int i=redshift_index; i<zstart_index; i++)
	{
		z = thisIntegralTable->zmin + i*thisIntegralTable->dz;
		tmp += integral_table[numz*numf*dcell_index + numz*factor_index + i];
// 		printf("z=%e:\t %e\n",z,tmp);
	}
	tmp = tmp*photHI*1.e12;
	
	return tmp;
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
		printf("Error: NO File!");
		exit(1);
	}
	fseek(fp, 0, SEEK_END);
	len_byte = ftell(fp);
	rewind(fp);
	
	num = len_byte/(sizeof(double)*3);
	data = malloc(num*3*sizeof(double));
	fread(data, sizeof(double), 3*num, fp);
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
	array = malloc(3*num_size*sizeof(double));
	
	for(int i=0; i<num_size; i++)
	{
		redshift = zmin + d*(i + offset_write);
		set_norm_pdf(pdf_params, redshift);
		array[3*i] = redshift;
		array[3*i+1] = pdf_params->amplitude;
		array[3*i+2] = pdf_params->constant;
	}
	
	write_table(num_size*3, offset_write*3, array, filename);
	
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

double *read_table_integrals(char *filename, integral_table_t *thisIntegralTable)
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
	fread(array, sizeof(double), num, fp);
	fclose(fp);
	
	numz = (thisIntegralTable->zmax - thisIntegralTable->zmin)/thisIntegralTable->dz+1;
	
	numf = (thisIntegralTable->fmax - thisIntegralTable->fmin)/thisIntegralTable->df+1;
	
	numdcell = (thisIntegralTable->dcellmax - thisIntegralTable->dcellmin)/thisIntegralTable->ddcell+1;
	
	assert(num == numz*numf*numdcell);
	
// 	for(int i=0; i<numdcell; i++) printf("%e\n",array[i]);
	
	return array;
}

integral_table_t * initIntegralTable(double zmin, double zmax, double dz, double fmin, double fmax, double df, double dcellmin, double dcellmax, double ddcell)
{
	integral_table_t *newIntegralTable;
	newIntegralTable = malloc(sizeof(integral_table_t));
	
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
	params = malloc(sizeof(dens_integrand_t));
	
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
	params = malloc(sizeof(dens_integrand_t));
	
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
	
	write_table(num1*num2*3, offset_write*num2*3, array, filename);
	
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
	
// 	printf("tmp = %e\t%e\t%e\t%e\n", tmp, dens_cell, x, 2.*SQR(frac*7.61/((1.+x)*dens_cell)));
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
	params = malloc(sizeof(redshift_integrand_t));
	
	params->simParam = simParam;
	
	for(int i=0; i<num1; i++)
	{
		dens_cell = pow(10.,thisRedshiftTable->dens_cell_min+thisRedshiftTable->ddens_cell*i);
		params->dens_cell = dens_cell;
		for(int j=0; j<num2; j++)
		{
			z = thisRedshiftTable->zmin+j*dz;
// 			printf("z = %e\t dens_cell = %e\n",z,dens_cell);
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
	params = malloc(sizeof(redshift_integrand_t));
	
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
	
	write_table(num1*num2*3, offset_write*num2*3, array, filename);
	
	free(params);
	free(array);
}

//------------------------------------------------------------------------------
// write table
//------------------------------------------------------------------------------

void write_table(int num, int offset, double *array, char *filename)
{
#ifdef __MPI
	MPI_File mpifile;
	MPI_Offset offset_mpi;
	MPI_Status status;
	
	offset_mpi = offset*sizeof(double);
	
	printf("num = %d\t offset = %d\n", num, offset);
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifile);
	MPI_File_write_at_all(mpifile, offset_mpi, array, num, MPI_DOUBLE, &status);
	MPI_File_close(&mpifile);
#else
	FILE * fp;
	
	fp = fopen(filename, "wb");
	fwrite(array, sizeof(double), num, fp);
	fclose(fp);
#endif
}