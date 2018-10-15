 #include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <dirent.h>

#ifdef __MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "utils.h"

int file_exist(char *filename)
{
	FILE *file;
	if((file = fopen (filename, "rt"))){
		fclose(file);
		return 1;
	}else{
		return 0;
	}
}

int directory_exist(char *dirname)
{
        DIR* dir = opendir(dirname);
        if(dir)
        {
                closedir(dir);
                return 1;
        }
        else
        {
                return 0;
        }
}

void *get_directory(char *filename)
{
        char *token = NULL;
        char *directory = NULL;
        size_t length;

        token = strrchr(filename, '/');

        if( filename == NULL )
        {
                printf("There is no file:'%s'\n", filename); /* You decide here */
        }

        if (token == NULL) 
        {
                printf("No directory indicated\n"); /* You decide here */
        }
        else
        {
            length = strlen(token);
            directory = malloc(length);
            memcpy(directory, token+1, length);
        }

        return directory;
}

// void abort()
// {
// #ifdef __MPI
//         MPI_Bcast(&error, 1, MPI_INT, );
//         MPI_Abort(MPI_COMM_WORLD);
//         
// #else
//         exit(EXIT_FAILURE);
// #endif
// }