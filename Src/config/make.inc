#
#   definition of the variables for the Makefile that makes the pw library.
#

#   If you are lucky the defaults for each type of CPU may work and give you a
#   fast functionning code.  But as that depends on the compilers and libraries
#   installed on your machine, nothing is assured.
#
#   You will find below alternatives that may work by switching commented
#   with uncommented lines.

#   Indicate what type of cpu you are using (uncomment one).  It will only suggest the best
#   compiler for the job.  "cat /proc/cpuinfo" will tell you what CPU is on your machine.

CPU=Intel
# CPU=AMD
# CPU=ARM


#   Suggestions for compilers.  Feel free to override the suggestions.
#   As of early 2024 ifx should be considered "beta" and unreliable.

ifeq ($(CPU),Intel)
   F77=ifort
   # F77=ifx
   # F77=gfortran
   # F77=pgfortran
else
   ifeq ($(CPU),AMD)
      F77=gfortran
      # F77=ifx
      # F77=ifort
      # F77=pgfortran
   else
      F77=gfortran
   endif
endif

#   debugging option (true, false).  Use if the code throws a segmentation fault and you
#   do not feel like migrating to Quantum Espresso...
#   Uncomment also optional changes for debugging below for fine tuning...

DEBUG=false
# DEBUG=true



#   below this point you should have a good reason (fine tune debugging, curiosity) to change
#   the defaults.  Each compiler has its own section.


#   Where esdf modules can be found

ESDFDIR=LIBS/esdf_read


#   suffix for executables if you want to compare compilers

SUFFIX:=_$(F77)
# SUFFIX=

ifeq ($(DEBUG),true)
   SUFFIX:=$(SUFFIX)_debug
endif


#   compatibility with older compilers, before fortran2008 (true/false)

PRE2008=false
# PRE2008=true



#  For Intel ifort and ifx compiler:

IFORTX=false
ifeq ($(F77),ifort)
   IFORTX=true
endif
ifeq ($(F77),ifx)
  IFORTX=true
endif

ifeq ($(IFORTX),true)

   #   profiling option (-g,empty).  Use it for tunning the code.

   PROF=
   # PROF=-g

   #   The next three choices indicate how FFTS and linear algebra are dealt during execution.

   #   indicates if blas and lapack must be compiled (true/false), instead of relying on libraries.

   COMPBLAS=false
   # COMPBLAS=true

   # it is safer to use fortran dcg...

   # MKLEXTRALIB=false
   MKLEXTRALIB=true

   #   Choice of FFT package (ncar/vfft/fftw3/mkldfti).

   # FFT=vfft
   # FFT=ncar
   # FFT=fftw3
   FFT=mkldfti

   #   Choice of hk_psi implementation (reference, mkldfti)

   HKPSI=mkldfti
   # HKPSI=reference

   ifeq ($(DEBUG),true)
      HKPSI=reference
   endif

   #   Compiling options and libraries.  If compiling for another machine beware of "-xHost"
   #   For an old fortran compiler use -mkl instead of -qmkl

   # LIBS=-mkl
   LIBS=-qmkl

   ifeq ($(DEBUG),false)
      FCOPTS=-xHost -O3 $(PROF)
      FCOMP=-qopenmp
      # FCOMP=-openmp
   else
      FCOPTS=-check all -warn all,nodec,interfaces,noexternal -gen-interfaces -traceback -fpe0 -fp-stack-check -g -O0
      FCOMP=
   endif

   FCEXTRA=
   FCSAFE=-O0

endif




#  For gcc gfortran compiler:


ifeq ($(F77),gfortran)

   #   profiling option (-pg,empty).  Use it for tunning the code.

   PROF=
   # PROF=-pg

   #   The next three choices indicate how FFTS and linear algebra are dealt during execution.

   #   indicates if blas and lapack must be compiled (true/false), instead of relying on libraries.

   COMPBLAS=false
   # COMPBLAS=true

   # MKLEXTRALIB=false
   MKLEXTRALIB=true

   ifeq ($(DEBUG),true)
      COMPBLAS=true
   endif

   #   Choice of FFT package (ncar/vfft/fftw3/mkldfti).
   #   Beware of recent implementations of fftw3.   Seems they changed APIs...

   FFT=vfft
   # FFT=ncar
   # FFT=fftw3
   # FFT=mkldfti

   #   Choice of hk_psi implementation (reference, mkldfti)

   # HKPSI=mkldfti
   HKPSI=reference

   ifeq ($(DEBUG),true)
      HKPSI=reference
   endif

   #   Compiling options and libraries.  Must have the chosen library installed.

   ifeq ($(DEBUG),false)

      FCOPTS=-fomit-frame-pointer -O3 -mtune=native $(PROF)
      FCOMP=-fopenmp
      FCEXTRA=-funroll-all-loops -ffast-math

      # choose one set of libraries.  openblas seems to be the best

      # defaults
      # LIBS=-lblas -llapack -lgomp

      # MKL
      # LIBS=-lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -lgomp

      # atlas  adjust LIBDIR
      # LIBDIR =/usr/local/lib
      # LIBS =-L$LIBDIR  -llapack -lf77blas -lcblas -latlas

      # openblas
      LIBS=-lopenblas -lpthread -lgomp

      # AMD-blis  adjust LIBDIR,  fftw3 has problems...
      # LIBDIR =/opt/AMD/aocl/aocl-linux-gcc-2.2.0/lib
      # LIBS=-L$LIBDIR -lflame -lblis
      # LIBS=-L$LIBDIR -lflame -lblis -lfftw3

   else

      FCOPTS=-g -Wall -fbounds-check -pedantic-errors -ffpe-trap=zero,invalid,overflow,underflow
      FCOMP=
      FCEXTRA=

      LIBS=

   endif

   FCSAFE=-O0

endif




#  For portland group pgfortran compiler, it is much less used than the previous, so beware.


ifeq ($(F77),pgfortran)

   #   Never tried profiling or debugging

   PROF=

   #   The next three choices indicate how FFTS and linear algebra are dealt during execution.

   #   indicates if blas and lapack must be compiled (true/false), instead of relying on libraries.

   COMPBLAS=false
   # COMPBLAS=true

   # MKLEXTRALIB=false
   MKLEXTRALIB=true

   #   Choice of FFT package (ncar/vfft/fftw3/mkldfti).
   #   Beware of recent implementations of fftw3.   Seems they changed APIs...

   FFT=vfft
   # FFT=ncar
   # FFT=fftw3

   #   Choice of hk_psi implementation (reference, mkldfti)

   HKPSI=reference

   #   Compiling options and libraries.  Must have the chosen library installed.

   FCOPTS=-fast
   FCOMP=
   FCEXTRA=

   # choose one set of libraries.  openblas seems to be the best.  Check also if flexiblas is available

   # defaults
   LIBS=-lblas -llapack

   # Maybe could be default
   # LIBS=-lflexiblas

   # AMD-blis  adjust LIBDIR,  fftw3 has problems...
   # LIBDIR =/opt/AMD/aocl/aocl-linux-gcc-2.2.0/lib
   # LIBS=-L$LIBDIR -lflame -lblis
   # LIBS=-L$LIBDIR -lflame -lblis -lfftw3


   FCSAFE=-O0

endif


