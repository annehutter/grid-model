#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "utils.h"
#include "confObj.h"
#include "sources.h" 

#include "nion_sources.h"

int main (int argc, /*const*/ char * argv[]) { 

    char iniFile[MAXLENGTH];
    confObj_t simParam;
    
    sourcelist_t *sourcelist = NULL;
    sourcelist_t *sourcelist_HeI = NULL;
    sourcelist_t *sourcelist_HeII = NULL;

    //parse command line arguments and be nice to user
    if (argc != 2) {
        printf("cifog: (C)  - Use at own risk...\n");
        printf("USAGE:\n");
        printf("cifog iniFile\n");
        
        exit(EXIT_FAILURE);
    } else {
        strcpy(iniFile, argv[1]);
    }
    
    //read paramter file
    simParam = readConfObj(iniFile);
    
    if(file_exist(simParam->sources_file) == 1){
        sourcelist = read_sources(simParam->sources_file);
        sourcelist_HeI = read_sources(simParam->sources_file);
        sourcelist_HeII = read_sources(simParam->sources_file);
    }else{
        fprintf(stderr, "Could not read source file %s\n", simParam->sources_file);
        exit(EXIT_FAILURE);
    }
    
    printf("sourcelist->numSources = %d\n", sourcelist->numSources);
    printf("sourcelist_HeI->numSources = %d\n", sourcelist_HeI->numSources);
    printf("sourcelist_HeII->numSources = %d\n", sourcelist_HeII->numSources);

    modify_nion_fesc(sourcelist, simParam->factorNion, simParam->factorFesc);
    compute_HeII_nion(simParam, sourcelist, sourcelist_HeI);
    compute_HeIII_nion(simParam, sourcelist, sourcelist_HeII);

    write_sources(sourcelist, simParam->sources_file_new);
    write_sources(sourcelist_HeI, simParam->sources_file_HeI);
    write_sources(sourcelist_HeII, simParam->sources_file_HeII);
    
    deallocate_sourcelist(sourcelist);
    deallocate_sourcelist(sourcelist_HeI);
    deallocate_sourcelist(sourcelist_HeII);

    confObj_del(&simParam);
}
