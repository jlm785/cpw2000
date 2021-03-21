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
!>     to a file. 
!>     It is commented for gnuplot or xmgrace!

       subroutine dos_print_file(filename,nhist,ehist,dhist,chist,       &
     &    lidos,ezero,efguess,vcell,ztot)

!      written December 7 , 2013 from old code. jlm
!      Modified, documentation, 19 September 2020. JLM
!      copyright  J.L.Martins, INESC-MN.

!      version 4.53 of cpw

       implicit none

       integer, parameter  :: REAL64 = selected_real_kind(12)

!      input:

       character(len=*) , intent(in)   ::  filename                      !<  input file name
       integer, intent(in)             ::  nhist                         !<  number of points in histograms
       real(REAL64), intent(in)        ::  ehist(nhist)                  !<  energies of the histogram
       real(REAL64), intent(in)        ::  dhist(nhist)                  !<  density of states
       real(REAL64), intent(in)        ::  chist(nhist)                  !<  integrated density of states
       logical, intent(in)             ::  lidos                         !<  true if integrated density of states is to be computed.
       real(REAL64), intent(in)        ::  ezero                         !<  position of internal neutral level
       real(REAL64), intent(in)        ::  efguess                       !<  estimate of the Fermi energy
       real(REAL64), intent(in)        ::  vcell                         !<  unit cell volume (in atomic units)
       real(REAL64), intent(in)        ::  ztot                          !<  total charge density (electrons/cell)

!      counters

       integer     :: i

!      constants

       real(REAL64), parameter :: HARTREE = 27.21138386_REAL64

!      printout to filename

       open(unit=15,file=filename,form='formatted')

       write(15,'("#")')
       write(15,'("# ",f14.5,"  ezero (localization of average ",        &
     &    "internal potential)")') -ezero
       write(15,'("# ",f14.5,"  efguess (Fermi level estimate)")')       &
     &          efguess
       write(15,'("# ",f14.5,"  vcell (cell volume a.u. )")') vcell
       write(15,'("# ",f14.5,"  ztot (electrons/cell )")') ztot
       write(15,'("# ",l5,"  Integrated DOS available")') lidos
       write(15,'("#")')
       write(15,'("#   Energy in eV")')
       write(15,'("#   Density of States  (DOS) in el/cell/eV")')
       if(lidos) then
         write(15,'("#   Integrated DOS in el/cell")')
         write(15,'("#")')
         do i=1,nhist
           write(15,'(2x,3f15.8)') ehist(i)*HARTREE,dhist(i)/HARTREE,    &
     &             chist(i)
         enddo
       else
         write(15,'("#")')
         do i=1,nhist
           write(15,'(2x,3f15.8)') ehist(i)*HARTREE,dhist(i)/HARTREE
         enddo
       endif

       close(unit=15)

       return
       end subroutine dos_print_file
