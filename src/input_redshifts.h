#ifndef INPUT_REDSHIFTS_H
#define INPUT_REDSHIFTS_H

double *read_redshift_list(char *redshift_file, int num_snapshots);
double *initRedshift_list(int num_snapshots);
void deallocateRedshift_list(double *redshift_list);

#endif