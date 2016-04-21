#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <complex.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "confObj.h"
#include "grid.h"

#define SQR(X) ((X) * (X))


void construct_tophat_filter(fftw_complex *filter, int nbins, ptrdiff_t local_0_start, ptrdiff_t local_n0, float smooth_scale)
{
	double normCoeff;
	const double half_nbins = nbins*0.5;
	const double sq_smooth_scale = 0.25*SQR(smooth_scale);
	
	normCoeff = 0;
	for(int i=0; i<nbins; i++){	  
		const double i_expr = half_nbins-abs(i - half_nbins);		
		const double sq_i_expr = SQR(i_expr);
		for(int j=0; j<nbins; j++){
			const double j_expr = half_nbins - abs(j - half_nbins);
			const double sq_j_expr = SQR(j_expr);
			for(int k=0; k<nbins; k++){
			  const double k_expr = half_nbins - abs(k - half_nbins);
			  const double sq_k_expr = SQR(k_expr);
			  
			  const double expr = sq_i_expr + sq_j_expr + sq_k_expr - sq_smooth_scale;
			  normCoeff += (expr<=0.) ? 1.:0.;
			}
		}
	}
	normCoeff = 1./normCoeff;
	
	for(int i=0; i<local_n0; i++)
	{
		const double i_expr = half_nbins-abs(i + local_0_start - half_nbins);		
		const double sq_i_expr = SQR(i_expr);
		for(int j=0; j<nbins; j++)
		{
		  	const double j_expr = half_nbins - abs(j - half_nbins);
			const double sq_j_expr = SQR(j_expr);
			for(int k=0; k<nbins; k++)
			{			
				const double k_expr = half_nbins - abs(k - half_nbins);
				const double sq_k_expr = SQR(k_expr);
			  
				const double expr = sq_i_expr + sq_j_expr + sq_k_expr - sq_smooth_scale;
				
				filter[i*nbins*nbins+j*nbins+k] = (expr<=0.) ? normCoeff+0.*I : 0.+0.*I;
			}
		}
	}
}

void convolve_fft(grid_t *thisGrid, fftw_complex *filter, fftw_complex *nion_smooth)
{
	int nbins;
	double factor;
#ifdef __MPI
	ptrdiff_t alloc_local, local_n0, local_0_start;
#else
	ptrdiff_t local_n0;
#endif
	
	fftw_complex *nion;
	fftw_complex *nion_ft, *filter_ft;
	fftw_plan plan_nion, plan_filter, plan_back;
	
	nbins = thisGrid->nbins;
	nion = thisGrid->frac_Q;
	local_n0 = thisGrid->local_n0;
	
#ifdef __MPI
	local_0_start = thisGrid->local_0_start;
	
	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
	assert(local_0_start == thisGrid->local_0_start);
	assert(local_n0 == thisGrid->local_n0);
	
	nion_ft = fftw_alloc_complex(alloc_local);
	filter_ft = fftw_alloc_complex(alloc_local);
	
	plan_nion = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, nion, nion_ft, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MPI_TRANSPOSED_OUT); //FFTW_MPI_TRANSPOSED_OUT
	plan_filter = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, filter, filter_ft, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MPI_TRANSPOSED_OUT);
#else 
	nion_ft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	filter_ft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	
	plan_nion = fftw_plan_dft_3d(nbins, nbins, nbins, nion, nion_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_filter = fftw_plan_dft_3d(nbins, nbins, nbins, filter, filter_ft, FFTW_FORWARD, FFTW_ESTIMATE);
#endif
	
	fftw_execute(plan_nion);
	fftw_execute(plan_filter);
	
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
	plan_back = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, nion_ft, nion_smooth, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MPI_TRANSPOSED_IN);
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
			}
		}
	}	
	
	fftw_destroy_plan(plan_nion);
	fftw_destroy_plan(plan_filter);
	fftw_destroy_plan(plan_back);
	
	fftw_free(nion_ft);
	fftw_free(filter_ft);
}

void determine_ion_fractions(fftw_complex *nion_smooth, int nbins, ptrdiff_t local_n0, int smallest_scale)
{
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				if(creal(nion_smooth[i*nbins*nbins+j*nbins+k])>=1.)
				{
					nion_smooth[i*nbins*nbins+j*nbins+k] = 1.+0.*I;
				}else{
					if(smallest_scale==0)
					{
						nion_smooth[i*nbins*nbins+j*nbins+k] = 0.+0.*I;
					}
				}
			}
		}
	}
}

void choose_ion_fraction(fftw_complex *nion_smooth, fftw_complex *XHII_tmp, grid_t *thisGrid)
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
				if(creal(nion_smooth[i*nbins*nbins+j*nbins+k])>=creal(XHII_tmp[i*nbins*nbins+j*nbins+k]))
				{
					XHII_tmp[i*nbins*nbins+j*nbins+k] =  nion_smooth[i*nbins*nbins+j*nbins+k];
				}
			}
		}
	}
}


void map_bubbles_to_nrec(fftw_complex *XHII_tmp, grid_t *thisGrid)
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
				thisGrid->nrec[i*nbins*nbins+j*nbins+k] = thisGrid->nrec[i*nbins*nbins+j*nbins+k]*XHII_tmp[i*nbins*nbins+j*nbins+k];
			}
		}
	}
}

void combine_bubble_and_web_model(fftw_complex *XHII_tmp, grid_t *thisGrid)
{
  	int nbins;
	ptrdiff_t local_n0;
	fftw_complex *XHII;
	
	nbins = thisGrid->nbins;
	local_n0 = thisGrid->local_n0;
	XHII = thisGrid->XHII;
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				XHII[i*nbins*nbins+j*nbins+k] = XHII[i*nbins*nbins+j*nbins+k]*XHII_tmp[i*nbins*nbins+j*nbins+k];
			}
		}
	}
}

void copy_grid_array(fftw_complex *XHII_tmp, grid_t *thisGrid)
{
  	int nbins;
	ptrdiff_t local_n0;
	fftw_complex *XHII;
	
	nbins = thisGrid->nbins;
	local_n0 = thisGrid->local_n0;
	XHII = thisGrid->XHII;
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				XHII[i*nbins*nbins+j*nbins+k] = XHII_tmp[i*nbins*nbins+j*nbins+k];
			}
		}
	}
}

void compute_ionization_field(confObj_t simParam, grid_t *thisGrid)
{
	int nbins;
	float smooth_scale;
	float box_size;
	float lin_scales, inc_log_scales;
	
	fftw_complex *filter;
	fftw_complex *nion_smooth;
	fftw_complex *XHII_tmp;
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
	
#ifdef __MPI
	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
	filter = fftw_alloc_complex(alloc_local);
	nion_smooth = fftw_alloc_complex(alloc_local);
	XHII_tmp = fftw_alloc_complex(alloc_local);
#else
	local_0_start = thisGrid->local_0_start;
	local_n0 = thisGrid->local_n0;
	filter = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	nion_smooth = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	XHII_tmp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
#endif
	initialize_grid(XHII_tmp, nbins, local_n0, 0.);
	
	
	lin_bins = lin_scales/box_size*(float)nbins;
	factor_exponent = inc_log_scales/lin_bins;
	num_scales = (log(box_size/lin_scales) + lin_bins*log(1.+factor_exponent))/log(1. + factor_exponent );
	printf("lin_bins = %e\t factor_exponent = %e\t num_scales = %d\n", lin_bins, factor_exponent, num_scales);
	
	for(int scale=0; scale<num_scales; scale++)
	{
	  	for(int i=0; i<local_n0; i++)
		{
			for(int j=0; j<nbins; j++)
			{
				for(int k=0; k<nbins; k++)
				{
					nion_smooth[i*nbins*nbins+j*nbins+k] = thisGrid->frac_Q[i*nbins*nbins+j*nbins+k];
				}
			}
		}
		
		
		inc = (float)num_scales - (float)scale;
		
		if(inc <= lin_bins) smooth_scale = inc;
		else smooth_scale = lin_bins*pow(1. + factor_exponent, inc-lin_bins);
		
		printf("inc = %e\t lin_bins = %e\t scale = %d\t smooth_scale = %e\n",inc, lin_bins, scale, smooth_scale);

		construct_tophat_filter(filter, nbins, local_0_start, local_n0, smooth_scale);
		
		convolve_fft(thisGrid, filter, nion_smooth);

		if(scale==num_scales-1)
		{
			determine_ion_fractions(nion_smooth, nbins, local_n0, 1);
		}else{
			determine_ion_fractions(nion_smooth, nbins, local_n0, 0);
		}
		
		choose_ion_fraction(nion_smooth, XHII_tmp, thisGrid);
	}
	
	if(simParam->use_web_model == 1)
	{
		combine_bubble_and_web_model(XHII_tmp, thisGrid);
		map_bubbles_to_nrec(XHII_tmp, thisGrid);
	}else{
		copy_grid_array(XHII_tmp, thisGrid);
	}

	fftw_free(filter);
	fftw_free(nion_smooth);
	fftw_free(XHII_tmp);
}


