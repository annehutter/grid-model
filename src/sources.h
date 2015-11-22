#ifndef SOURCES_H
#define SOURCES_H
#endif

typedef struct
{
	float pos[3];
	double Nion;
	float fesc;
} source_t;


/*functions*/
source_t *source_init();
void source_free(source_t *thisSource);
int count_sources(char *filename);
source_t *allocate_source_list(int num_sources);
void deallocate_source_list(source_t thisSource[]);
void read_sources(char * filename, int num_sources, source_t *thisSource);