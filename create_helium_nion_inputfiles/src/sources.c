#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "sources.h"
#include "utils.h"

/* functions on sourcelist struct (list of sources) ---------------------------------------------------------------*/

sourcelist_t *allocate_sourcelist(int Nallocated)
{
	sourcelist_t *newSourcelist;
	newSourcelist = malloc(sizeof(sourcelist_t));
	if(newSourcelist == NULL)
	{
		fprintf(stderr, "newSourcelist in sourcelist_init (sources.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	newSourcelist->numSources = 0;
	newSourcelist->Nallocated = Nallocated;
	newSourcelist->source = malloc(Nallocated*sizeof(source_t));
	if(newSourcelist->source == NULL)
	{
		fprintf(stderr, "newSourcelist->source in sourcelist_init (sources.c) could not be allocated\n");
		exit(EXIT_FAILURE);
	}
	
	return newSourcelist;
}

void deallocate_sourcelist(sourcelist_t *thisSourcelist)
{	
	if(thisSourcelist != NULL)
	{
		free(thisSourcelist->source);
		free(thisSourcelist);
	}
}

sourcelist_t *read_sources(char * filename)
{
	float xpos, ypos, zpos;
	float fesc;
	double Nion;
	int spectraID;
	
	int numSources;

	FILE * fp;
	char line[MAXLENGTH];
	
	fp = fopen(filename, "rt");
	if(fp==NULL)
	{
		fprintf(stderr, "Can not open file %s\n",filename);
		exit(EXIT_FAILURE);
	}
	
	if(fgets(line,MAXLENGTH,fp)!=NULL)
	{
		sscanf(line,"%d\n",&numSources);
	}
        
	sourcelist_t *newSourcelist = allocate_sourcelist(numSources);

	while(fgets(line,MAXLENGTH,fp)!=NULL)
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
	assert(numSources == newSourcelist->numSources);
	numSources = newSourcelist->numSources;
	newSourcelist->source = realloc(newSourcelist->source, numSources*sizeof(source_t));
	
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

void write_sources(sourcelist_t *thisSourcelist, char *filename)
{
    int numSources = thisSourcelist->numSources;
    source_t *thisSource;

    FILE *fp;
    fp = fopen(filename, "w");
    if(fp==NULL)
	{
		fprintf(stderr, "Can not open file %s\n",filename);
		exit(EXIT_FAILURE);
	}
	
    fprintf(fp, "%d\n", numSources);

	for(int i=0; i<numSources; i++)
    {
        thisSource = &(thisSourcelist->source[i]);
        
        fprintf(fp, "%f\t%f\t%f\t%e\t%d\t%f\n", thisSource->pos[0], thisSource->pos[1], thisSource->pos[2], thisSource->Nion, i, thisSource->fesc);
    }
	
	fclose(fp);
}
 
