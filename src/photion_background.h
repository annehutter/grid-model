#ifndef PHOTION_BACKGROUND_H
#define PHOTION_BACKGROUND_H

typedef struct
{
    int num;
    int Nallocated;
    double *redshift;
    double *photHI;
    double *QHII;
} photIonlist_t;


photIonlist_t *allocate_photIonlist(int Nallocated);
void deallocate_photIonlist(photIonlist_t *thisPhotIonlist);
photIonlist_t *read_photIonlist(char *filename);
double get_photHI_from_redshift(const photIonlist_t *thisPhotIonlist, double redshift);
double get_photHI_from_fillingfactor(const photIonlist_t *thisPhotIonlist, double QHII);

#endif