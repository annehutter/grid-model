#ifndef RECOMBINATION
#define RECOMBINATION
#endif

typedef struct
{
	double dens;
	double photHI;
	double temp;
	double redshift;
	confObj_t simParam;
	
	double amplitude;
	double constant;
	double beta;
} nrec_params_t;

void set_nrec_params(nrec_params_t *nrec_params, pdf_params_t *pdf_params, confObj_t simParam, double dens, double photHI, double temp, double redshift);
double nrec_rate(double x, void * p);
double calc_nrec_rate(nrec_params_t *params, double redshift);
double calc_nrec_rate_integrand(double x, void * p);
double calc_nrec_integral(nrec_params_t *params, double zstart, double zend);
double get_nrec_history(nrec_params_t *nrec_params, double zstart, double redshift);
double get_nrec_timestep(nrec_params_t *nrec_params, double zstart, double zend);
void compute_number_recobinations(grid_t *thisGrid, confObj_t simParam);
