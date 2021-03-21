!------------------------------------------------------------!
! This file is distributed as part of the cpw2000 code and   !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the cpw2000        !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the cpw2000 code is not yet written         !
!                                                            !
! The cpw2000 code is hosted on GitHub:                      !
!                                                            !
! https://github.com/jlm785/cpw2000                          !
!------------------------------------------------------------!

!>     writes the density of states and integrated density of states 
!>     to a file to be plotted by gnuplot.

       subroutine dos_out_gnuplot(ioreplay,filename,                     &
     &               nhist,ehist,dhist,chist,lidos,lask)

!      adapted June 12 , 2014. JLM
!      Modified, documentation, 19 September 2020. JLM
!      copyright  J.L.Martins, INESC-MN.

!      version 4.53 of cpw

       implicit none

       integer, parameter  :: REAL64 = selected_real_kind(12)

!      input:

       integer, intent(in)             :: ioreplay                       !<  tape number for reproducing calculations

       character(len=*) , intent(in)   ::  filename                      !<  input file
       integer, intent(in)             ::  nhist                         !<  number of points in histograms
       real(REAL64), intent(in)        ::  ehist(nhist)                  !<  energies of the histogram
       real(REAL64), intent(in)        ::  dhist(nhist)                  !<  density of states
       real(REAL64), intent(in)        ::  chist(nhist)                  !<  integrated density of states
       logical, intent(in)             ::  lidos                         !<  true if integrated density of states is to be computed.
       logical, intent(in)             ::  lask                          !<  if true asks interactively for a plot to be displayed

!      local variables

       real(REAL64)      ::  ar, br
       real(REAL64)      ::  emin,emax
       character(len=1)  ::  yesno
       integer           ::  ixe, io

!      counters

       integer     ::  i

!      constants

       real(REAL64), parameter :: HARTREE = 27.21138386_REAL64


!      finds nice plotrange

       ar = (ehist(nhist)-ehist(1))*HARTREE / 10

       call plot_step(ar,br)

       br = br/HARTREE

       ixe = int(ehist(1)/br)
       emin = br*(ixe-1)*HARTREE
       ixe = int(ehist(nhist)/br)
       emax = br*(ixe)*HARTREE

!      printout to filename

       io = 15
       open(unit=io,file=filename,form='formatted')


       write(io,*) "set terminal wxt enhanced"
       write(io,*)
       write(io,*) "set title  'Density of States' font 'Helvetica-Bold' "
       write(io,*) "set ylabel 'DOS (el/cell/eV)' font  'Helvetica-Bold' "
       write(io,*) "set xlabel 'Energy (eV)' font 'Helvetica-Bold' "
       write(io,*) "set xrange [ ",emin," : ",emax," ] " 
       write(io,*) "set border lw 3"
       write(io,*) "set xtics font 'Helvetica-Bold' "
       write(io,*) "set ytics font 'Helvetica-Bold' "
       write(io,*) "plot '-'  u 1:2 w lines title 'DOS'; pause -1 "
       write(io,*)
       if(lidos) then
         do i=1,nhist
           write(io,'(2x,3f15.8)') ehist(i)*HARTREE,dhist(i)/HARTREE,    &
     &             chist(i)
         enddo
       else
         do i=1,nhist
           write(io,'(2x,3f15.8)') ehist(i)*HARTREE,dhist(i)/HARTREE
         enddo
       endif
       write(io,*)

       close(unit=io)

       write(6,*)
       write(6,*) '  The density of States plot is in ',filename
       write(6,*)

!      asks whether a figure should be displayed on screen

       if(lask) then
         write(6,*)
         write(6,*) '  Do you want to see the plot now? (y/n)'
         write(6,*)
       
         read(5,*) yesno
         write(ioreplay,*) yesno,'   see plot'

         if(yesno == 'y' .or. yesno == 'Y') then

           write(6,*)
           write(6,*) '  Hit return to continue '
           write(6,*)

!           call system("gnuplot " // filename // "  2> /dev/null ")
           call execute_command_line("gnuplot " // filename // "  2> /dev/null ")

         endif

       endif
      
       return
       end subroutine dos_out_gnuplot
