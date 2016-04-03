#ifndef MEAN_FREE_PATH
#define MEAN_FREE_PATH
#endif

typedef struct
{
	double redshift;
	double amplitude;
	double constant;
	double beta;
} pdf_params_t;

double Hubble(confObj_t simParam);
double pdf(double x, void * p);
double calc_integral(pdf_params_t params, double upLim);
double pdf_mass(double x, void * p);
double calc_integral_mass(pdf_params_t params, double upLim);
void set_norm_pdf(pdf_params_t * params, confObj_t simParam);
double frac_densSS(double densSS, confObj_t simParam);
double calc_mfp(confObj_t simParam);