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
		./src/cifog.c
OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := cifog

# USE-MPI = YES

OPTIMIZE = -O3 -ftree-vectorize
WARNING = -Wall -Wextra -Wshadow -Wpedantic -g

FFTW3DIR :=/opt/local/include
FFTW_CFLAGS := -I$(FFTW3DIR)
FFTW3LIBDIR :=/opt/local/lib
FFTW3_LINK := -L$(FFTW3LIBDIR) -lfftw3 -Xlinker -rpath -Xlinker $(FFTW3_LIBDIR) 

GSL_FOUND := $(shell gsl-config --version)
ifndef GSL_FOUND
  $(error $(ccred)Error:$(ccreset) GSL not found in path - please install GSL before installing $(DISTNAME).$(VERSION) $(ccreset))
endif
GSL_CFLAGS := $(shell gsl-config --cflags)
GSL_LIBDIR := $(shell gsl-config --prefix)/lib
GSL_LINK := $(shell gsl-config --libs) -Xlinker -rpath -Xlinker $(GSL_LIBDIR)

LDFLAGS := $(GSL_LINK) $(FFTW3_LINK) -lm
CFLAGS := -c -std=c99 -march=native $(WARNING) $(OPTIMIZE) $(GSL_CFLAGS) $(FFTW_CFLAGS)

UNAME := $(shell uname) 

ifeq ($(UNAME), Darwin)
	COMPILER := clang
else
	COMPILER := gcc
endif

ifdef USE-MPI
	CC := mpicc
	CFLAGS += -D __MPI
	LDFLAGS += -lmpich
        FFTW3_LINK +=  -lfftw3_mpi
else
	CC := $(COMPILER)
endif

.PHONY: all clean clena celan celna

all: $(SOURCES) $(EXECUTABLE)

celan celna clena:clean


$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -rf $(OBJECTS) $(DOBJECTS) $(EXECUTABLE)
