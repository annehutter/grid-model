#include <stdio.h>
#include <stdlib.h>

#include "sources.h"

/* functions on source struct (for one source) ---------------------------------------------------------------*/

source_t *source_init()
{
	source_t *newSource;
	newSource = malloc(sizeof(source_t));
	if(newSource == NULL)
	{
		fprintf(stderr, "ERROR: source_init Not enough memory to allocate source. \n");
		exit(EXIT_FAILURE);
	}
	
	for(int i=0; i<3; i++) newSource->pos[i] = 0.;
	newSource->Nion = 0.;
	newSource->fesc = 0.;
}

void source_free(source_t *thisSource)
{
	free(thisSource);
}

/* functions on sourcee data ------------------------------------------------------------------*/

int count_sources(char *filename)
{
	int number=0;
	int ch=0;
	
	FILE * fp;
	fp = fopen(filename, "r");
	if(fp==NULL)
	{
		fprintf(stderr, "Can not open file %s\n",filename);
		exit(EXIT_FAILURE);
	}

	while(!feof(fp))
	{
		ch = fgetc(fp);
		if(ch == '\n') number++;
	}
	
	fclose(fp);

	return number;
}

source_t *allocate_source_list(int num_sources)
{
	source_t *newSource;
	newSource = malloc(sizeof(source_t)*num_sources);
	if(newSource == NULL)
	{
		fprintf(stderr, "ERROR: source_init Not enought memory to allocate particle.\n");
		exit(EXIT_FAILURE);
	}
	
	for(int source=0; source<num_sources; source++)
	{
		newSource[source].Nion = 0.;
		newSource[source].fesc = 0.;
		for(int i=0; i<3; i++) newSource[source].pos[i] = 0.;
	}
	
	return newSource;
}

void deallocate_source_list(source_t thisSource[])
{
	free(thisSource);
}

/* reading source list -----------------------------------------------------------------------------*/

void read_sources(char * filename, int num_sources, source_t *thisSource)
{
// 	source_t *thisSource;
// 	int num_sources;

	float xpos, ypos, zpos;
	float fesc;
	double Nion;
	int spectraID;
	
	int num;
	int source;
	
// 	num_sources = count_sources(filename);
// 	thisSource = allocate_source_list(num_sources);

	FILE * fp;
	char line[256];
	
	fp = fopen(filename, "rt");
	if(fp==NULL)
	{
		fprintf(stderr, "Can not open file %s\n",filename);
		exit(EXIT_FAILURE);
	}
	
	if(fgets(line,256,fp)!=NULL)
	{
		sscanf(line,"%d\n",&num);
	}

	source = 0;
	while(fgets(line,256,fp)!=NULL)
	{
		sscanf(line,"%f\t%f\t%f\t%lf\t%d\t%f\n", &xpos, &ypos, &zpos, &Nion, &spectraID, &fesc);
		thisSource[source].Nion = Nion;
		thisSource[source].fesc = fesc;
		thisSource[source].pos[0] = xpos;
		thisSource[source].pos[1] = ypos;
		thisSource[source].pos[2] = zpos;	
// 		printf("%e\t%f\n", thisSource[source].Nion,thisSource[source].fesc);
		source++;
	}
	  
	fclose(fp);
}