SOURCES := 	main.c \
		confObj.c \
		parse_ini.c \
		xmem.c \
		xstring.c \
		grid.c \
		sources.c \
		sources_to_grid.c \
		fraction_q.c \
		filtering.c \
		phys_const.c \
		self_shielding.c \
		density_distribution.c \
		recombination.c \
		mean_free_path.c \
		convolution_fftw.c \
		utils.c\
		input_redshifts.c\
		input_grid.c \
		photion_background.c \
		redshift_tools.c \
		cifog.c

SOURCES := $(addprefix src/, $(SOURCES))
OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := cifog

#USE-MPI = YES

include common.mk

.PHONY: all clean clena celan celna

all: $(SOURCES) $(EXECUTABLE)

celan celna clena:clean


$(EXECUTABLE): $(OBJECTS)
	$(COMPILER) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(COMPILER) $(CFLAGS) $< -o $@

clean: 
	rm -rf $(OBJECTS) $(DOBJECTS) $(EXECUTABLE)
