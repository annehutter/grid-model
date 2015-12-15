#ifndef GRID_H
#define GRID_H
#endif

/* structure for grid */

typedef struct
{
	int nbins;
	float box_size;
	
	float xmin, ymin, zmin;
	
	fftw_complex *igm_density;
	fftw_complex *halo_density;
	fftw_complex *igm_clump;
	fftw_complex *nion;
	
	fftw_complex *XHII;
	
	fftw_complex *photHI;
	double mean_photHI;
	
	int local_n0;
	int local_0_start;
} grid_t;

/* functions */

grid_t *initGrid();
void read_files_to_grid(grid_t *thisGrid, confObj_t thisInput);
void read_grid(fftw_complex *toThisArray, int nbins, int local_n0, int local_0_start, char *filename);
void initialize_grid(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, double value);
void deallocate_grid(grid_t *thisGrid);
void write_grid_to_file_float(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename);
void write_grid_to_file_double(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename);
void save_to_file_XHII(grid_t *thisGrid, char *filename);
void save_to_file_photHI(grid_t *thisGrid, char *filename);