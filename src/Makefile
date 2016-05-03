SOURCES := 	./main.c \
		./confObj.c \
		./parse_ini.c \
		./xmem.c \
		./xstring.c \
		./grid.c \
		./sources.c \
		./sources_to_grid.c \
		./fraction_q.c \
		./filtering.c \
		./phys_const.c \
		./self_shielding.c \
		./density_distribution.c \
		./recombination.c \
		./mean_free_path.c
OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := cifog

USE-MPI = YES

OPTIMIZE = -O3 -ftree-vectorize
WARNING = -Wall -Wextra -Wshadow -g

ifdef USE-MPI
	CC := mpicc
	CFLAGS := -c -std=c99 -march=native -lm $(WARNING) $(OPTIMIZE) -D __MPI
	LDFLAGS := -lfftw3_mpi -lfftw3 -lm -lmpich -lgsl -lgslcblas

else
	CC := gcc
	CFLAGS := -c -std=c99 -march=native -lm $(WARNING) $(OPTIMIZE) 
	LDFLAGS := -lfftw3 -lm -lgsl -lgslcblas
endif

.PHONY: all clean

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -rf $(OBJECTS) $(DOBJECTS) $(EXECUTABLE)
