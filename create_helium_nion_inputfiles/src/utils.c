 #include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <complex.h>

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
