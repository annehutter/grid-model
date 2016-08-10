#ifndef GRID_H
#define GRID_H
#endif

/* structure for grid */

typedef struct
{
	int nbins;
	float box_size;
	float lin_scales;
	float inc_log_scales;
	
	float xmin, ymin, zmin;
	
	fftw_complex *igm_density;
    fftw_complex *igm_clump;
	fftw_complex *nion;
	fftw_complex *cum_nion;
	fftw_complex *cum_nabs;
	fftw_complex *frac_Q;
	
	fftw_complex *XHII;
	fftw_complex *nrec;
	
	fftw_complex *photHI;
	double mean_photHI;
	
	int local_n0;
	int local_0_start;
} grid_t;

/* functions */

grid_t *initGrid();
void read_files_to_grid(grid_t *thisGrid, confObj_t thisInput);
void read_nion(grid_t *thisGrid, char *filename, int double_precision);
void read_igm_density(grid_t *thisGrid, char *filename, int double_precision);
void read_igm_clump(grid_t *thisGrid, char *filename, int double_precision);

#ifdef __MPI
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename);
void read_grid_doubleprecision(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename);
#else
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, char *filename);
void read_grid_doubleprecision(fftw_complex *toThisArray, int nbins, int local_n0, char *filename);
#endif

void initialize_grid(fftw_complex *thisArray, int nbins, int local_n0, double value);
void deallocate_grid(grid_t *thisGrid);

#ifdef __MPI
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename);
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename);
#else
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, char *filename);
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, char *filename);
#endif

void save_to_file_XHII(grid_t *thisGrid, char *filename);
void save_to_file_photHI(grid_t *thisGrid, char *filename);
