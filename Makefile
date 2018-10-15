SOURCES := 	./src/main.c \
		./src/confObj.c \
		./src/parse_ini.c \
		./src/xmem.c \
		./src/xstring.c \
		./src/grid.c \
		./src/sources.c \
		./src/sources_to_grid.c \
		./src/fraction_q.c \
		./src/filtering.c \
		./src/phys_const.c \
		./src/self_shielding.c \
		./src/density_distribution.c \
		./src/recombination.c \
		./src/mean_free_path.c \
		./src/convolution_fftw.c \
		./src/utils.c\
		./src/input_redshifts.c\
		./src/input_grid.c \
		./src/photion_background.c \
		./src/redshift_tools.c \
		./src/cifog.c \
		./src/checks.c

OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := cifog

USE-MPI=YES
 
include common.mk


.PHONY: all clean clena celan celna

all: $(SOURCES) $(EXECUTABLE)

celan celna clena:clean


$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -rf $(OBJECTS) $(DOBJECTS) $(EXECUTABLE)
