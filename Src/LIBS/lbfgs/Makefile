ROOT = ../..

include $(ROOT)/make.inc

src77 = dpofa.f  dtrsl.f  lbfgs.f  setulb.f  timer.f



obj =  $(src77:.f=.o)

libadd:    $(obj)  dpmeps.o
	ar -rv $(ROOT)/libpw$(SUFFIX).a *.o

dpmeps.o: dpmeps.f
	$(F77) -c $(FCSAFE) dpmeps.f

%.o: %.f
	$(F77) -c $(FCOPTS)  $<

clean:
	rm -f *.o *.mod

.PHONY:   clean libadd
