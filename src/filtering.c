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

#include "fraction_q.h"
#include "convolution_fftw.h"

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

// versatile function (thisGrid is only used for grid dimensions and domain decomposition
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

// versatile function (thisGrid is only used for grid dimensions and domain decomposition
void choose_ion_fraction(fftw_complex *nion_smooth, fftw_complex *Xion_tmp, grid_t *thisGrid)
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
				if(creal(nion_smooth[i*nbins*nbins+j*nbins+k])>=creal(Xion_tmp[i*nbins*nbins+j*nbins+k]))
				{
					Xion_tmp[i*nbins*nbins+j*nbins+k] =  nion_smooth[i*nbins*nbins+j*nbins+k];
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

void compute_ionization_field(confObj_t simParam, grid_t *thisGrid, int specie)
{
	int nbins;
	float smooth_scale;
	float box_size;
	float lin_scales, inc_log_scales;
	    
	fftw_complex *filter;
	fftw_complex *nion_smooth;
    fftw_complex *nabs_smooth;
    fftw_complex *frac_Q_smooth;
	fftw_complex *Xion_tmp;
    
    fftw_complex *cum_nion;
    fftw_complex *cum_nabs;
    fftw_complex *nrec;
    fftw_complex *Xion;
    
    if(specie == 1){
        cum_nion = thisGrid->cum_nion_HeI;
        cum_nabs = thisGrid->cum_nabs_HeI;
        nrec = thisGrid->nrec_HeI;
        Xion = thisGrid->XHeII;
    }else if(specie == 2){
        cum_nion = thisGrid->cum_nion_HeII;
        cum_nabs = thisGrid->cum_nabs_HeII;
        nrec = thisGrid->nrec_HeII;
        Xion = thisGrid->XHeIII;
    }else{
        cum_nion = thisGrid->cum_nion;
        cum_nabs = thisGrid->cum_nabs;
        nrec = thisGrid->nrec;
        Xion = thisGrid->XHII;
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
	
#ifdef __MPI
	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
	filter = fftw_alloc_complex(alloc_local);
	nion_smooth = fftw_alloc_complex(alloc_local);
    nabs_smooth = fftw_alloc_complex(alloc_local);
    frac_Q_smooth = fftw_alloc_complex(alloc_local);
    
	Xion_tmp = fftw_alloc_complex(alloc_local);
#else
	local_0_start = thisGrid->local_0_start;
	local_n0 = thisGrid->local_n0;
	filter = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	if(filter == NULL)
	{
		fprintf(stderr, "filter in compute_ionization_field (filtering.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	nion_smooth = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	if(nion_smooth == NULL)
	{
		fprintf(stderr, "nion_smooth in compute_ionization_field (filtering.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	nabs_smooth = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
	if(nabs_smooth == NULL)
	{
		fprintf(stderr, "nabs_smooth in compute_ionization_field (filtering.c) could not be allocated\n");
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
#endif
	initialize_grid(Xion_tmp, nbins, local_n0, 0.);
	
	
	lin_bins = lin_scales/box_size*(float)nbins;
	factor_exponent = inc_log_scales/lin_bins;
	num_scales = (log(box_size/lin_scales) + lin_bins*log(1.+factor_exponent))/log(1. + factor_exponent );
	printf("\n #linear bins = %e\t factor_exponent = %e\t #tophat filter sizes = %d\n", lin_bins, factor_exponent, num_scales);
	
	for(int scale=0; scale<num_scales; scale++)
	{
	  	for(int i=0; i<local_n0; i++)
		{
			for(int j=0; j<nbins; j++)
			{
				for(int k=0; k<nbins; k++)
				{
					nion_smooth[i*nbins*nbins+j*nbins+k] = creal(cum_nion[i*nbins*nbins+j*nbins+k])+0.*I;
                    nabs_smooth[i*nbins*nbins+j*nbins+k] = creal(cum_nabs[i*nbins*nbins+j*nbins+k])+0.*I;
				}
			}
		}
		
		
		inc = (float)num_scales - (float)scale;
		
		if(inc <= lin_bins) smooth_scale = inc;
		else smooth_scale = lin_bins*pow(1. + factor_exponent, inc-lin_bins);
		
		printf("  inc = %e\t lin_bins = %e\t scale = %d\t smooth_scale = %e bins\n",inc, lin_bins, scale, smooth_scale);

		construct_tophat_filter(filter, nbins, local_0_start, local_n0, smooth_scale);
		
		convolve_fft(thisGrid, filter, nion_smooth, cum_nion);
		convolve_fft(thisGrid, filter, nabs_smooth, cum_nabs);

        compute_Q(thisGrid, frac_Q_smooth, nion_smooth, nabs_smooth);
        
		if(scale==num_scales-1)
		{
			determine_ion_fractions(frac_Q_smooth, nbins, local_n0, 1);
		}else{
			determine_ion_fractions(frac_Q_smooth, nbins, local_n0, 0);
		}
		
		choose_ion_fraction(frac_Q_smooth, Xion_tmp, thisGrid);
	}
	
	
	if(simParam->use_web_model == 1 && (specie != 1 && specie !=2))
	{
        printf("\n\ncombining web and bubbles!\n\n");
		combine_bubble_and_web_model(Xion_tmp, Xion, thisGrid);
		map_bubbles_to_nrec(Xion_tmp, nrec, thisGrid);
	}else{
		copy_grid_array(Xion_tmp, Xion, thisGrid);
	}

	fftw_free(filter);
	fftw_free(nion_smooth);
    fftw_free(nabs_smooth);
    fftw_free(frac_Q_smooth);
	fftw_free(Xion_tmp);
}

 
// void compute_ionization_field(confObj_t simParam, grid_t *thisGrid)
// {
// 	int nbins;
// 	float smooth_scale;
// 	float box_size;
// 	float lin_scales, inc_log_scales;
// 	
// 	fftw_complex *filter;
// 	fftw_complex *nion_smooth;
// 	fftw_complex *XHII_tmp;
// #ifdef __MPI
// 	ptrdiff_t alloc_local, local_n0, local_0_start;
// #else
// 	ptrdiff_t local_n0, local_0_start;
// #endif
// 	
// 	float inc, factor_exponent, lin_bins;
// 	int num_scales;
// 
// 	nbins = thisGrid->nbins;
// 	box_size = thisGrid->box_size;
// 	lin_scales = thisGrid->lin_scales;
// 	inc_log_scales = thisGrid->inc_log_scales;
// 	
// #ifdef __MPI
// 	alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
// 	filter = fftw_alloc_complex(alloc_local);
// 	nion_smooth = fftw_alloc_complex(alloc_local);
// 	XHII_tmp = fftw_alloc_complex(alloc_local);
// #else
// 	local_0_start = thisGrid->local_0_start;
// 	local_n0 = thisGrid->local_n0;
// 	filter = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
// 	if(filter == NULL)
// 	{
// 		fprintf(stderr, "filter in compute_ionization_field (filtering.c) could not be allocated\n");
// 		exit(EXIT_FAILURE);
// 	}
// 	nion_smooth = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
// 	if(nion_smooth == NULL)
// 	{
// 		fprintf(stderr, "nion_smooth in compute_ionization_field (filtering.c) could not be allocated\n");
// 		exit(EXIT_FAILURE);
// 	}
// 	XHII_tmp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
// 	if(XHII_tmp == NULL)
// 	{
// 		fprintf(stderr, "XHII_tmp in compute_ionization_field (filtering.c) could not be allocated\n");
// 		exit(EXIT_FAILURE);
// 	}
// #endif
// 	initialize_grid(XHII_tmp, nbins, local_n0, 0.);
// 	
// 	
// 	lin_bins = lin_scales/box_size*(float)nbins;
// 	factor_exponent = inc_log_scales/lin_bins;
// 	num_scales = (log(box_size/lin_scales) + lin_bins*log(1.+factor_exponent))/log(1. + factor_exponent );
// 	printf("\n #linear bins = %e\t factor_exponent = %e\t #tophat filter sizes = %d\n", lin_bins, factor_exponent, num_scales);
// 	
// 	for(int scale=0; scale<num_scales; scale++)
// 	{
// 	  	for(int i=0; i<local_n0; i++)
// 		{
// 			for(int j=0; j<nbins; j++)
// 			{
// 				for(int k=0; k<nbins; k++)
// 				{
// 					nion_smooth[i*nbins*nbins+j*nbins+k] = creal(thisGrid->frac_Q[i*nbins*nbins+j*nbins+k])+0.*I;
// 				}
// 			}
// 		}
// 		
// 		
// 		inc = (float)num_scales - (float)scale;
// 		
// 		if(inc <= lin_bins) smooth_scale = inc;
// 		else smooth_scale = lin_bins*pow(1. + factor_exponent, inc-lin_bins);
// 		
// 		printf("  inc = %e\t lin_bins = %e\t scale = %d\t smooth_scale = %e\n",inc, lin_bins, scale, smooth_scale);
// 
// 		construct_tophat_filter(filter, nbins, local_0_start, local_n0, smooth_scale);
// 		
// // 		convolve_fft_XHII(thisGrid, filter, nion_smooth);
// 		convolve_fft(thisGrid, filter, nion_smooth, thisGrid->frac_Q);
// 
// 		if(scale==num_scales-1)
// 		{
// 			determine_ion_fractions(nion_smooth, nbins, local_n0, 1);
// 		}else{
// 			determine_ion_fractions(nion_smooth, nbins, local_n0, 0);
// 		}
// 		
// 		choose_ion_fraction(nion_smooth, XHII_tmp, thisGrid);
// 	}
// 	
// 	if(simParam->use_web_model == 1)
// 	{
// 		combine_bubble_and_web_model(XHII_tmp, thisGrid);
// 		map_bubbles_to_nrec(XHII_tmp, thisGrid);
// 	}else{
// 		copy_grid_array(XHII_tmp, thisGrid);
// 	}
// 
// 	fftw_free(filter);
// 	fftw_free(nion_smooth);
// 	fftw_free(XHII_tmp);
// }
