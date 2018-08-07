OPTIMIZE = -O3 -ftree-vectorize
WARNING = -Wall -Wextra -Wshadow -g

FFTW3DIR :=/usr/local/include
FFTW_CFLAGS := -I$(FFTW3DIR)
FFTW3LIBDIR :=/usr/local/lib64
FFTW3_LINK := -L$(FFTW3LIBDIR) -lfftw3 -Xlinker -rpath -Xlinker $(FFTW3_LIBDIR) 

GSL_FOUND := $(shell gsl-config --version)
ifndef GSL_FOUND
  $(error $(ccred)Error:$(ccreset) GSL not found in path - please install GSL before installing $(DISTNAME).$(VERSION) $(ccreset))
endif
GSL_CFLAGS := $(shell gsl-config --cflags)
GSL_LIBDIR := $(shell gsl-config --prefix)/lib64
GSL_LINK := $(shell gsl-config --libs) -Xlinker -rpath -Xlinker $(GSL_LIBDIR)

LDFLAGS := -lm $(GSL_LINK) $(FFTW3_LINK)
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
    LDFLAGS += -lmpi -lfftw3_mpi
else
    CC := $(COMPILER)
endif
