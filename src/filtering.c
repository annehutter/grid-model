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
// 			  if(pow((float)(nbins/2-abs(i-nbins/2)),2)+pow((float)(nbins/2-abs(j-nbins/2)),2)+pow((float)(nbins/2-abs(k-nbins/2)),2)<=pow(smooth_scale*0.5,2)){
// 					normCoeff += 1.;
// 				}
			}
		}
	}
	normCoeff = 1./normCoeff;
	
// 	printf("normCoeff = %e\n",normCoeff);

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
// 				if(pow((float)(nbins/2-abs(i+local_0_start-nbins/2)),2)+pow((float)(nbins/2-abs(j-nbins/2)),2)+pow((float)(nbins/2-abs(k-nbins/2)),2)<=pow(smooth_scale*0.5,2))
// 				{
// 					filter[i*nbins*nbins+j*nbins+k] = normCoeff+0.*I;
// 				}else{
// 					filter[i*nbins*nbins+j*nbins+k] = 0.+0.*I;
// 				}
			}
		}
	}
}

void convolve_fft(grid_t *thisGrid, fftw_complex *filter, fftw_complex *nion_smooth)
{
	int nbins;
	double factor;
	ptrdiff_t alloc_local, local_n0, local_0_start;
	
// 	ptrdiff_t n0, n1, n2, local_n1, local_1_start, local_n0x, local_0x_start;
	
	fftw_complex *nion;
	fftw_complex *nion_ft, *filter_ft;
	fftw_plan plan_nion, plan_filter, plan_back;
	
	nbins = thisGrid->nbins;
	nion = thisGrid->nion;
	local_0_start = thisGrid->local_0_start;
	local_n0 = thisGrid->local_n0;
	
// 	n0 = nbins;
// 	n1 = nbins;
// 	n2 = nbins;
	
#ifdef __MPI
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
	
// 	fftw_mpi_local_size_3d_transposed(n0, n1, n2, MPI_COMM_WORLD, &local_n0x, &local_0x_start, &local_n1, &local_1_start);
	
// 	printf("local_n0 = %d\tlocal_0_start = %d\tlocal_n0x = %d\tlocal_0x_start = %d\tlocal_n1 = %d\tlocal_1_start = %d\n",local_n0,local_0_start,local_n0x,local_0x_start,local_n1,local_1_start);
	
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

void determine_ion_fractions(fftw_complex *nion_smooth, int nbins, ptrdiff_t local_0_start, ptrdiff_t local_n0, int smallest_scale)
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

void choose_ion_fraction(fftw_complex *nion_smooth, grid_t *thisGrid)
{
	int nbins;
	ptrdiff_t local_0_start, local_n0;
	fftw_complex *XHII;
	
	nbins = thisGrid->nbins;
	local_0_start = thisGrid->local_0_start;
	local_n0 = thisGrid->local_n0;
	XHII = thisGrid->XHII;
	
	for(int i=0; i<local_n0; i++)
	{
		for(int j=0; j<nbins; j++)
		{
			for(int k=0; k<nbins; k++)
			{
				if(creal(nion_smooth[i*nbins*nbins+j*nbins+k])>creal(XHII[i*nbins*nbins+j*nbins+k]))
				{
					XHII[i*nbins*nbins+j*nbins+k] = nion_smooth[i*nbins*nbins+j*nbins+k];
				}
			}
		}
	}
}

void compute_ionization_field(grid_t *thisGrid)
{
	int nbins;
	int half_nbins;
	float smooth_scale;
	fftw_complex *filter;
	fftw_complex *nion_smooth;
	ptrdiff_t alloc_local, local_n0, local_0_start;

	nbins = thisGrid->nbins;
	half_nbins = nbins*0.5;
	
#ifdef __MPI
	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
	filter = fftw_alloc_complex(alloc_local);
	nion_smooth = fftw_alloc_complex(alloc_local);
#else
	local_0_start = thisGrid->local_0_start;
	local_n0 = thisGrid->local_n0;
	filter = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	nion_smooth = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
#endif

	
	for(int scale=0; scale<nbins; scale++)
	{
		printf("scale = %d\n",scale);
		smooth_scale = nbins - (float)scale;
		construct_tophat_filter(filter, nbins, local_0_start, local_n0, smooth_scale);
		if(scale==126) write_grid_to_file_float(filter, thisGrid->nbins, thisGrid->local_n0, thisGrid->local_0_start, "filter.out");
		
		convolve_fft(thisGrid, filter, nion_smooth);
		if(scale==nbins-1)
		{
			determine_ion_fractions(nion_smooth, nbins, local_0_start, local_n0, 1);
		}else{
			determine_ion_fractions(nion_smooth, nbins, local_0_start, local_n0, 0);
		}
		
		choose_ion_fraction(nion_smooth, thisGrid);
	}
	
	fftw_free(filter);
	fftw_free(nion_smooth);
}


