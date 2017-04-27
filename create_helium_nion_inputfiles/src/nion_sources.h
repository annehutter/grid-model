#ifndef NION_SOURCES_H
#define NION_SOURCES_H
#endif

void modify_nion_fesc(sourcelist_t *thisSourcelist, double factorNion, double factorFesc);
void compute_HeII_nion(confObj_t simParam, sourcelist_t *thisSourcelist, sourcelist_t *thisSourcelist_HeI);
void compute_HeIII_nion(confObj_t simParam, sourcelist_t *thisSourcelist, sourcelist_t *thisSourcelist_HeIII);
