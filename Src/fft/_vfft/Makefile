ROOT = ../..

include $(ROOT)/make.inc

src77 = $(wildcard *.f)



obj =  $(src77:.f=.o)

libadd:    $(obj)
	ar -rv $(ROOT)/libpw$(SUFFIX).a *.o

%.o: %.f
	$(F77) -c $(FCOPTS) $<

clean:
	rm -f *.o *.mod

.PHONY:   clean libadd
