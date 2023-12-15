ROOT = ../..

include $(ROOT)/make.inc



libadd:   hk_psi_c16.o hk_psi_spin_c16.o
	ar -rv $(ROOT)/libpw$(SUFFIX).a *.o

hk_psi_c16.o: hk_psi_c16.f90
	$(F77) -c $(FCOPTS) $(FCEXTRA) $(FCOMP) hk_psi_c16.f90

hk_psi_spin_c16.o: hk_psi_spin_c16.f90
	$(F77) -c $(FCOPTS) $(FCEXTRA) $(FCOMP) hk_psi_spin_c16.f90

clean:
	rm -f *.o *.mod *genmod*

.PHONY:   clean libadd
