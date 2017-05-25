SOURCES := 	./src/main.c \
		./src/sources.c \
		./src/utils.c \
		./src/nion_sources.c \
		./src/confObj.c \
		./src/parse_ini.c \
		./src/xmem.c \
		./src/xstring.c \
		./src/phys_const.c \
		./src/cross_sections.c \
		./src/recombination_coeff.c
OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := create_helium_inputfiles

# USE-MPI = YES

OPTIMIZE = -O3 -ftree-vectorize
WARNING = -Wall -Wextra -Wshadow -g

ifdef USE-MPI
	CC := mpicc
	CFLAGS := -c -std=c99 -march=native -lm $(WARNING) $(OPTIMIZE)
	LDFLAGS := -lm 

else
	CC := gcc
	CFLAGS := -c -std=c99 -march=native -lm $(WARNING) $(OPTIMIZE) 
	#-D DEBUG_NREC
	LDFLAGS := -lm 
endif

.PHONY: all clean

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -rf $(OBJECTS) $(DOBJECTS) $(EXECUTABLE)
