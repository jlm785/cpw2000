
include make.inc

all: cpw.exe cpw_post_process.exe

cpw.exe:  libpw  cpw_variables.mod cpw.f90
	$(F77) -o cpw$(SUFFIX).exe $(FCOPTS) cpw.f90 libpw$(SUFFIX).a $(LIBS)

cpw_post_process.exe:  libpw  cpw_variables.mod cpw_post_process.f90
	$(F77) -o cpw_post_process$(SUFFIX).exe $(FCOPTS) cpw_post_process.f90 libpw$(SUFFIX).a $(LIBS)

libpw:     fftlib  timelib  hkpsilib    esdflib  blas   mklextra
	(cd cpw               &&  make libadd )
	(cd cpw_lib           &&  make libadd )
	(cd cpw_pp_lib        &&  make libadd )
	(cd LIBS/lbfgs        &&  make libadd )
	(cd fft/_cfft         &&  make libadd )
	(cd omp               &&  make libadd )
	(cd ao_interpolation  &&  make libadd )
	(cd ALFA              &&  make libadd )
	(cd crystal           &&  make libadd )
	(cd ditsp             &&  make libadd )
	(cd mesh              &&  make libadd )
	(cd numeric           &&  make libadd )
	(cd read_print        &&  make libadd )
	(cd QE                &&  make libadd )


cpw_variables.mod:  cpw_variables.f90
	$(F77) -c $(FCOPTS) cpw_variables.f90

blas:
	if [ $(COMPBLAS) = true ]; then \
	    (cd LIBS/blas     &&  make libadd ) ;\
	fi
	if [ $(COMPBLAS) = true ]; then \
	    (cd LIBS/lapack   &&  make libadd ) ;\
	fi

mklextra:
	if [ $(MKLEXTRALIB) = true ]; then \
	    (cd LIBS/mklextra &&  make libadd ) ;\
	fi

esdflib:
	(cd $(ESDFDIR)        &&  make libadd )

fftlib:
	if [ $(FFT) = ncar ]; then \
	   (cd fft/_ncar      &&  make libadd ) ;\
	fi
	if [ $(FFT) = ncar ]; then \
	   (cd LIBS/ncar      &&  make libadd ) ;\
	fi
	if [ $(FFT) = fftw3 ]; then \
	   (cd fft/_fftw3     &&  make libadd ) ;\
	fi
	if [ $(FFT) = mkldfti ]; then \
	   (cd fft/_mkl_dfti_2D  &&  make libadd ) ;\
	fi
	if [ $(FFT) = vfft ]; then \
	   (cd fft/_vfft      &&  make libadd ) ;\
	fi
	if [ $(FFT) = vfft ]; then \
	   (cd LIBS/vfft      &&  make libadd ) ;\
	fi

timelib:
	if [ $(PRE2008) = true ]; then \
	   (cd time/dummy     &&  make libadd );\
	fi
	if [ $(PRE2008) = false ]; then \
	   (cd time           &&  make libadd );\
	fi

hkpsilib:
	if [ $(HKPSI) = mkldfti ]; then \
	   (cd hk_psi/mkl_3D  &&  make libadd ) ;\
	fi
	if [ $(HKPSI) = reference ]; then \
	   (cd hk_psi/reference  &&  make libadd ) ;\
	fi

clean:
	(cd cpw               &&  make clean )
	(cd cpw_lib           &&  make clean )
	(cd cpw_pp_lib        &&  make clean )
	(cd LIBS/lbfgs        &&  make clean )
	(cd fft/_cfft         &&  make clean )
	(cd omp               &&  make clean )
	(cd ao_interpolation  &&  make clean )
	(cd ALFA              &&  make clean )
	(cd crystal           &&  make clean )
	(cd ditsp             &&  make clean )
	(cd mesh              &&  make clean )
	(cd numeric           &&  make clean )
	(cd read_print        &&  make clean )
	(cd QE                &&  make clean )
	(cd $(ESDFDIR)        &&  make clean )
	(cd LIBS/blas         &&  make clean )
	(cd LIBS/lapack       &&  make clean )
	(cd LIBS/ncar         &&  make clean )
	(cd LIBS/vfft         &&  make clean )
	(cd LIBS/mklextra     &&  make clean )
	(cd fft/_ncar         &&  make clean )
	(cd fft/_fftw3        &&  make clean )
	(cd fft/_mkl_dfti_2D  &&  make clean )
	(cd fft/_vfft         &&  make clean )
	(cd time              &&  make clean )
	(cd time/dummy        &&  make clean )
	(cd back_compat       &&  make clean )
	(cd hk_psi/mkl_3D     &&  make clean )
	(cd hk_psi/reference  &&  make clean )
	(rm *.mod)

.PHONY:   clean libpw  blas  esdflib  fftlib  timelib  hkpsilib  mklextra
