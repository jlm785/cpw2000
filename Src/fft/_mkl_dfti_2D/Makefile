ROOT = ../..

include $(ROOT)/make.inc

src90 = cfft_finish.f90  cfft_mlt_c16.f90  cfft_prepare.f90



obj =  $(src90:.f90=.o)

libadd:    $(obj) mkl_dfti.o  mkl_handle.o
	ar -rv $(ROOT)/libpw$(SUFFIX).a *.o

mkl_dfti.o:  mkl_dfti.f90
	$(F77) -c $(FCOPTS) mkl_dfti.f90

mkl_handle.o:  mkl_handle.f90 mkl_dfti.o
	$(F77) -c $(FCOPTS) mkl_handle.f90

%.o: %.f90 mkl_dfti.o mkl_handle.o
	$(F77) -c $(FCOPTS)  -I mkl_dfti.h $<

clean:
	rm -f *.o *.mod

.PHONY:   clean libadd
