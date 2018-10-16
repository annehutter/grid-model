#ifndef RESTART_H
#define RESTART_H

int32_t save_restart_file(confObj_t simParam, grid_t *grid, int32_t cycle, int32_t snap, int32_t myRank);
int32_t read_restart_file(confObj_t simParam, grid_t *grid, int32_t *start_cycle, int32_t *snap);

#endif