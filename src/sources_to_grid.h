#ifndef  SOURCES_TO_GRID_H
#define SOURCES_TO_GRID_H

void read_update_nion(confObj_t simParam, sourcelist_t *thisSourcelist, grid_t *thisGrid, int snap);
void read_update_nion_HeI(confObj_t simParam, sourcelist_t *thisSourcelist, grid_t *thisGrid, int snap);
void read_update_nion_HeII(confObj_t simParam, sourcelist_t *thisSourcelist, grid_t *thisGrid, int snap);

void map_nion_to_grid(fftw_complex *thisNionArray, grid_t *thisGrid, sourcelist_t *thisSourcelist);

#endif