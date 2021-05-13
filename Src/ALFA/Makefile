ROOT = ..

include $(ROOT)/make.inc

src90 = $(wildcard *.f90)



obj =  $(src90:.f90=.o)

libadd: mods   $(obj)
	ar -rv $(ROOT)/libpw$(SUFFIX).a *.o

mods:
	$(F77) -c $(FCOPTS) $(ROOT)/cpw_variables.f90
	$(F77) -c $(FCOPTS) $(ROOT)/$(ESDFDIR)/esdf_key_mod.f90
	$(F77) -c $(FCOPTS) $(ROOT)/$(ESDFDIR)/esdf_mod.f90

%.o: %.f90
	$(F77) -c $(FCOPTS) $<

clean:
	rm -f *.o *.mod *genmod*

.PHONY:   clean libadd mods
