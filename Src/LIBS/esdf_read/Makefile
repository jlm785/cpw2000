ROOT = ../..

include $(ROOT)/make.inc



libadd:    esdf_key_mod.o  esdf_mod.o
	 ar -rv $(ROOT)/libpw$(SUFFIX).a *.o

esdf_key_mod.o: esdf_key_mod.f90
	$(F77) -c $(FCOPTS) esdf_key_mod.f90

esdf_mod.o: esdf_mod.f90   esdf_key_mod.o
	$(F77) -c $(FCOPTS) esdf_mod.f90

clean:
	rm -f *.o *.mod

.PHONY:   clean libadd
