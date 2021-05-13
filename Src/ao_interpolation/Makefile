ROOT = ..

include $(ROOT)/make.inc

src90 = ao_h_and_s.f90                           ao_interpolation_sub.f90          \
        ao_h_and_s_spin_orbit.f90                ao_interpolation_w90.f90          \
        ao_interpolation_core.f90                ao_interpolation_write_chk.f90    \
        ao_interpolation_out_band.f90            ao_interpolation_write_hr.f90     \
        ao_interpolation_out_band_fold_full.f90  ao_interpolation_write_wt_in.f90  \
        ao_interpolation_out_dos.f90             fourier_interpolation.f90         \
        ao_interpolation_out_ie.f90              irredbzinit.f90                   \
        ao_interpolation_prepare.f90             nonorthointerp.f90



obj =  $(src90:.f90=.o)

libadd:    $(obj)
	ar -rv $(ROOT)/libpw$(SUFFIX).a *.o

mods: $(ROOT)/cpw_variables.f90
	$(F77) -c $(FCOPTS) $(ROOT)/cpw_variables.f90

fourier_interpolation.o:  fourier_interpolation.f90
	$(F77) -c $(FCOPTS) fourier_interpolation.f90

nonorthointerp.o:  nonorthointerp.f90  fourier_interpolation.o
	$(F77) -c $(FCOPTS) nonorthointerp.f90

irredbzinit.o:  irredbzinit.f90
	$(F77) -c $(FCOPTS) irredbzinit.f90

ao_interpolation_sub.o:  ao_interpolation_sub.f90 mods nonorthointerp.o
	$(F77) -c $(FCOPTS) ao_interpolation_sub.f90

ao_h_and_s.o:  ao_h_and_s.f90
	$(F77) -c $(FCOPTS) ao_h_and_s.f90

ao_h_and_s_spin_orbit.o:  ao_h_and_s_spin_orbit.f90
	$(F77) -c $(FCOPTS) ao_h_and_s_spin_orbit.f90

ao_interpolation_core.o:  ao_interpolation_core.f90
	$(F77) -c $(FCOPTS) ao_interpolation_core.f90

ao_interpolation_out_band.o:  ao_interpolation_out_band.f90 nonorthointerp.o
	$(F77) -c $(FCOPTS) ao_interpolation_out_band.f90

ao_interpolation_out_band_fold_full.o:  ao_interpolation_out_band_fold_full.f90 nonorthointerp.o
	$(F77) -c $(FCOPTS) ao_interpolation_out_band_fold_full.f90

ao_interpolation_out_dos.o:  ao_interpolation_out_dos.f90 nonorthointerp.o
	$(F77) -c $(FCOPTS) ao_interpolation_out_dos.f90

ao_interpolation_out_ie.o:  ao_interpolation_out_ie.f90 nonorthointerp.o
	$(F77) -c $(FCOPTS) ao_interpolation_out_ie.f90

ao_interpolation_prepare.o:  ao_interpolation_prepare.f90 nonorthointerp.o irredbzinit.o
	$(F77) -c $(FCOPTS) ao_interpolation_prepare.f90

ao_interpolation_w90.o:  ao_interpolation_w90.f90 nonorthointerp.o
	$(F77) -c $(FCOPTS) ao_interpolation_w90.f90

ao_interpolation_write_chk.o:  ao_interpolation_write_chk.f90 nonorthointerp.o
	$(F77) -c $(FCOPTS) ao_interpolation_write_chk.f90

ao_interpolation_write_hr.o:  ao_interpolation_write_hr.f90 nonorthointerp.o
	$(F77) -c $(FCOPTS) ao_interpolation_write_hr.f90

ao_interpolation_write_wt_in.o:  ao_interpolation_write_wt_in.f90 nonorthointerp.o
	$(F77) -c $(FCOPTS) ao_interpolation_write_wt_in.f90


clean:
	rm -f *.o *.mod *genmod*

.PHONY:   clean libadd mods
