#ifndef CIFOG_H
#define CIFOG_H
#endif

int cifog(confObj_t simParam, const double *redshift_list, grid_t *grid, sourcelist_t *sourcelist,
          const integral_table_t *integralTable, photIonlist_t *photIonBgList, const int num_cycles, const int myRank);
