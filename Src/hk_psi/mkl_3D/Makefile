ROOT = ../..

include $(ROOT)/make.inc



libadd:   mkl_dfti.o  hk_psi_c16.o
	ar -rv $(ROOT)/libpw$(SUFFIX).a *.o

mkl_dfti.o: mkl_dfti.f90
	$(F77) -c $(FCOPTS)  mkl_dfti.f90

hk_psi_c16.o: hk_psi_c16.f90   mkl_dfti.o
	$(F77) -c $(FCOPTS) $(FCEXTRA) $(FCOMP) hk_psi_c16.f90

clean:
	rm -f *.o *.mod

.PHONY:   clean libadd
