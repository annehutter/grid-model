#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "sources.h"

// /* functions on source struct (for one source) ---------------------------------------------------------------*/
// 
// source_t *source_init()
// {
// 	source_t *newSource;
// 	newSource = malloc(sizeof(source_t));
// 	if(newSource == NULL)
// 	{
// 		fprintf(stderr, "ERROR: source_init Not enough memory to allocate source. \n");
// 		exit(EXIT_FAILURE);
// 	}
// 	
// 	for(int i=0; i<3; i++) newSource->pos[i] = 0.;
// 	newSource->Nion = 0.;
// 	newSource->fesc = 0.;
// }
// 
// void deallocate_source(source_t *thisSource)
// {
// 	assert(thisSource != NULL);
// 	free(thisSource);
// }

/* functions on sourcelist struct (list of sources) ---------------------------------------------------------------*/

sourcelist_t *allocate_sourcelist(int Nallocated)
{
	sourcelist_t *newSourcelist;
	newSourcelist = malloc(sizeof(sourcelist_t));
	if(newSourcelist == NULL)
	{
		fprintf(stderr, "ERROR: sourcelist_init: Not enough memory.\n");
		exit(EXIT_FAILURE);
	}
	
	newSourcelist->numSources = 0;
	newSourcelist->Nallocated = Nallocated;
	newSourcelist->source = malloc(Nallocated*sizeof(source_t));
	if(newSourcelist->source == NULL)
	{
		fprintf(stderr, "ERROR: sourcelist->source_init: Not enough memory.\n");
		exit(EXIT_FAILURE);
	}
	
// 	for(int i=0; i<numSources; i++)
// 	{
// 		source_t *newSource;
// 		newSource = malloc(sizeof(source_t));
// 		if(newSource == NULL)
// 		{
// 			fprintf(stderr, "ERROR: source_init: Not enough memory.\n");
// 			exit(EXIT_FAILURE);
// 		}
// 		newSource->Nion = 0.;
// 		newSource->fesc = 0.;
// 		for(int i=0; i<3; i++) newSource->pos[i] = 0.;
// 		newSourcelist->source[i] = newSource;
// 	}
	
	return newSourcelist;
}

// void deallocate_source_list(source_t thisSource[])
// {
// 	free(thisSource);
// }

void deallocate_sourcelist(sourcelist_t *thisSourcelist)
{	
	assert(thisSourcelist != NULL);
	
// 	numSources = thisSourcelist->numSources;
// 	
// 	for(int i=0; i<numSources; i++)
// 	{
// 		deallocate_source(thisSourcelist->source[i]);
// 	}
	
	free(thisSourcelist->source);
	free(thisSourcelist);
}

sourcelist_t *read_sources(char * filename)
{
	float xpos, ypos, zpos;
	float fesc;
	double Nion;
	int spectraID;
	
	int numSources;

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
		sscanf(line,"%d\n",&numSources);
	}
		
	sourcelist_t *newSourcelist = allocate_sourcelist(numSources);

	while(fgets(line,256,fp)!=NULL)
	{
		sscanf(line,"%f\t%f\t%f\t%lf\t%d\t%f\n", &xpos, &ypos, &zpos, &Nion, &spectraID, &fesc);
		
		source_t *newSource = &(newSourcelist->source[newSourcelist->numSources]);
		
		newSource->Nion = Nion;
		newSource->fesc = fesc;
		newSource->pos[0] = xpos;
		newSource->pos[1] = ypos;
		newSource->pos[2] = zpos;
		
		newSourcelist->numSources++;
		
		if(newSourcelist->numSources == newSourcelist->Nallocated)
		{
		  int nallocated = newSourcelist->Nallocated*1.1;
		  while(nallocated == newSourcelist->Nallocated) nallocated++;
		    
		  newSourcelist->source = realloc(newSourcelist->source, nallocated*sizeof(source_t));
		  if(newSourcelist->source == NULL)
		  {
			  fprintf(stderr, "ERROR: additional sources not allocatable\n");
			  exit(EXIT_FAILURE);
		  }
			newSourcelist->Nallocated = nallocated;
		}
	}
	numSources = newSourcelist->numSources;
	newSourcelist->source = realloc(newSourcelist->source, numSources*sizeof(source_t*));
	
	fclose(fp);
	
	return newSourcelist;
}

/* functions on source data ------------------------------------------------------------------*/

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


 
