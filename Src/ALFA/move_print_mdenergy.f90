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

!>     prints the md related energies

       subroutine move_print_mdenergy(flgcal,istep,                      &
     &      energy,ekin,ekcell,elects)

!      written 18 September 2002. jlm
!      modified (elects) 2 October 2013. jlm
!      Modified, documentation, June 2019. JLM
!      Modified, EPILNG, August 2019. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94 of cpw
!      version 1.5 of md


       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input:
       character(len=6), intent(in)    ::  flgcal                       !<  type of calculation
       integer, intent(in)             ::  istep                        !<  molecular dynamics step
       real(REAL64), intent(in)        ::  energy                       !<  kohn-sham energy in Hartree (including ewald)
       real(REAL64), intent(in)        ::  ekin                         !<  kinetic energy in Hartree of atoms (ions) 
       real(REAL64), intent(in)        ::  ekcell                       !<  kinetic energy in Hartree of cell (fictitious) 
       real(REAL64), intent(in)        ::  elects                       !<  electronic temperature*entropy in Hartree

       if(flgcal == 'MICRO ' .or. flgcal == 'LANG  ') then

         write(6,*)
         write(6,*)
         write(6,'(3x,i5,3f15.6,"  iteration, kohn-sham,",              &
     &        " atomickinetic, and totalenergy")')                      &
     &                  istep,energy,ekin,energy+ekin
         write(6,*)
         if(abs(elects) > 0.0001) write(6,'(3x,i5,f15.6,3x,             &
     &     "  iteration and free_energy")') istep,energy+ekin-elects
         write(6,*)

       elseif(flgcal == 'VCSLNG' .or. flgcal == 'VCSMIC' .or.           &
     &        flgcal == 'EPILNG') then

         write(6,*)
         write(6,*)
         write(6,'(3x,i6,3x,4f15.6,"  iteration, kohn-sham,",           &
     &        " atomickinetic, cellkinetic, and totalenergy")')         &
     &                istep,energy,ekin,ekcell,energy+ekin+ekcell
         write(6,*)
         if(abs(elects) > 0.0001) write(6,'(3x,i5,f15.6,3x,             &
     &     "  iteration and free_energy")')                             &
     &                   istep,energy+ekin+ekcell-elects
         write(6,*)

       else

         write(6,*)
         write(6,*)
         write(6,'(3x,i6,3x,f15.6,"  iteration, kohnsham energy",       &
     &        " = totalenergy")') istep,energy
         write(6,*)
         if(abs(elects) > 0.0001) write(6,'(3x,i5,f15.6,3x,             &
     &     "  iteration and free_energy")') istep,energy-elects
         write(6,*)

       endif

       return

       end subroutine move_print_mdenergy
