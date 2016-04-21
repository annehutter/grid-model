#ifndef SOURCES_H
#define SOURCES_H
#endif

typedef struct
{
	double Nion;
	float pos[3];
	float fesc;
} source_t;

typedef struct
{
	int numSources;
	int Nallocated;
	source_t *source;
} sourcelist_t;

/*functions*/
// source_t *source_init();
// void deallocate_source(source_t *thisSource);
sourcelist_t *allocate_sourcelist(int num_sources);
void deallocate_sourcelist(sourcelist_t *thisSourcelist);
sourcelist_t *read_sources(char * filename);
int count_sources(char *filename); 
