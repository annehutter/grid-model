CC = mpicc
CFLAGS = -c -std=c99 -D __MPI -march=native -O3
LDFLAGS = -lfftw3_mpi -lfftw3 -lm -lmpi 
SOURCES = 	./main.c \
		./confObj.c \
		./parse_ini.c \
		./xmem.c \
		./xstring.c \
		./grid.c \
		./sources.c \
		./sources_to_grid.c \
		./fraction_q.c \
		./filtering.c \
		./phys_const.c
OBJECTS = $(SOURCES:.c=.o)
DOBJECTS = $(SOURCES:.c=.d)
EXECUTABLE = cifog

.PHONY: all clean

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -rf $(OBJECTS) $(DOBJECTS) $(EXECUTABLE)
