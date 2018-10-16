#ifndef CIFOG_H
#define CIFOG_H

int cifog_init(char *iniFile, confObj_t *simParam, double **redshift_list, grid_t **grid,  integral_table_t **integralTable, photIonlist_t **photIonBgList, int *num_cycles, const int restart, const int myRank);

int cifog(confObj_t simParam, const double *redshift_list, grid_t *grid, sourcelist_t *sourcelist, const integral_table_t *integralTable, photIonlist_t *photIonBgList, const int num_cycles, const int myRank);
int cifog_step(confObj_t simParam, grid_t *grid, sourcelist_t *sourcelist, const integral_table_t *integralTable, photIonlist_t *photIonBgList, const int cycle, int snap, const int myRank);

int cifog_deallocate(confObj_t simParam, double *redshift_list, grid_t *grid, integral_table_t *integralTable, photIonlist_t *photIonBgList, const int myRank);
#endif