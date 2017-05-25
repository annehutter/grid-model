### Set the default compiler -- possible options are icc/gcc/clang
COMPILER:=

#### Add any compiler specific flags you want
CFLAGS:=

#### Add any compiler specific link flags you want
LDFLAGS:=


UNAME := $(shell uname)
## Colored text output
## Taken from: http://stackoverflow.com/questions/24144440/color-highlighting-of-makefile-warnings-and-errors
## Except, you have to use "echo -e" on linux and "echo" on Mac
ifeq ($(UNAME), Darwin)
  ECHO_COMMAND := echo
else
  ECHO_COMMAND := echo -e
endif
ifeq ($(TRAVIS_OS_NAME), linux)
  ECHO_COMMAND := echo
endif 
ccred:=$(shell $(ECHO_COMMAND) "\033[0;31m")
ccmagenta:=$(shell $(ECHO_COMMAND) "\033[0;35m")
ccgreen:=$(shell $(ECHO_COMMAND) "\033[0;32m")
ccblue:=$(shell $(ECHO_COMMAND) "\033[0;34m")
ccreset:=$(shell $(ECHO_COMMAND) "\033[0;0m")
boldfont:=$(shell $(ECHO_COMMAND) "\033[1m")
## end of colored text output


DO_CHECKS := 1
CLEAN_CMDS := celan celna clean clena distclean realclean
ifneq ($(filter $(CLEAN_CMDS),$(MAKECMDGOALS)),)
  DO_CHECKS := 0
endif


ifeq ($(DO_CHECKS), 1)
  ifeq ($(COMPILER),) 
    ## Set the default compiler (if not set)
    ifdef USE-MPI
      COMPILER := mpicc
    else
      ifeq ($(UNAME), Darwin)
        COMPILER := clang
      else
        COMPILER := gcc
      endif
    endif
  endif

  ## Check the hostname
  HOSTNAME := $(shell hostname)
  G2_STRING := hpc.swin.edu.au
  G2_FFTW_PKG := fftw/x86_64/gnu/3.3.3-openmpi-psm
  ON_G2 := 0
  ifeq ($(findstring $(G2_STRING),$(HOSTNAME)), $(G2_STRING))
    ON_G2 := 1
  endif

  ## Set the FFTW3 package
  ifneq (icc,$(findstring icc,$(COMPILER)))
    FFTW3DIR :=/opt/local/include
    FFTW3_CFLAGS := -I$(FFTW3DIR)
    FFTW3_LIBDIR :=/opt/local/lib
    FFTW3_LINK := -L$(FFTW3_LIBDIR) -lfftw3 -Xlinker -rpath -Xlinker $(FFTW3_LIBDIR)

    ifeq ($(ON_G2), 1)
      $(info On $(ccgreen)g2$(ccreset). Using package $(ccblue)$(G2_FFTW_PKG)$(ccreset))
      FFTW3DIR := /usr/local/x86_64/gnu/fftw-3.3.3-openmpi-psm
      FFTW3_CFLAGS := -I$(FFTW3DIR)
      FFTW3_LIBDIR := $(FFTW3DIR)/lib
      FFTW3_LINK := -L$(FFTW3_LIBDIR) -lfftw3 -Xlinker -rpath -Xlinker $(FFTW3_LIBDIR)
    endif
  else
    ## Using the intel compiler. Use the FFTW3 library from MKL
    $(info Using the intel compiler)
    FFTW3DIR := $(MKLROOT)
    FFTW3_LIBDIR := $(FFTW3DIR)/lib
    FFTW3_CFLAGS := -DMKL_ILP64 -I$(MKLROOT)/include/fftw
    FFTW3_LINK := -L$(FFTW3DIR)/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -Xlinker -rpath -Xlinker $(MKLROOT)/lib/intel64
  endif

  ## Use the fftw3 mpi library
  ifdef USE-MPI
    ifeq (icc,$(findstring icc,$(COMPILER)))
      FFTW3_LINK += -lmkl_blacs_intelmpi_ilp64
    else
      FFTW3_LINK += -lfftw3_mpi
    endif
  endif
  ## Done setting the FFTW3 package

  ## Set the GSL package
  GSL_FOUND := $(shell gsl-config --version)
  ifndef GSL_FOUND
    $(error $(ccred)Error:$(ccreset) GSL not found in path - please install GSL before installing $(DISTNAME).$(VERSION) $(ccreset))
  endif
  GSL_CFLAGS := $(shell gsl-config --cflags)
  GSL_LIBDIR := $(shell gsl-config --prefix)/lib
  GSL_LINK := $(shell gsl-config --libs) -Xlinker -rpath -Xlinker $(GSL_LIBDIR)
  ## Done setting the GSL library


  OPTIMIZE := -O3 -march=native
  WARNING := -Wall -Wextra -Wshadow -g
  LDFLAGS += $(GSL_LINK) $(FFTW3_LINK)
  ## Add the math library for other compilers 
  ifeq (icc,$(findstring icc,$(COMPILER)))
    OPTIMIZE += -xhost -ipo
  else
    OPTIMIZE += -ftree-vectorize -flto
    LDFLAGS += -lm
  endif
  CFLAGS := -c -std=gnu99 $(WARNING) $(OPTIMIZE) $(GSL_CFLAGS) $(FFTW3_CFLAGS)
  ifdef USE-MPI
    CFLAGS += -D__MPI
    LDFLAGS += $(FFTW_MPI_LINK)
  endif
endif ## end of DO_CHECKS -> do not run anything if make clean is running
