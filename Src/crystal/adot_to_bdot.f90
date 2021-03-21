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

!>     computes volume and reciprocal space metric

       subroutine adot_to_bdot(adot,vcell,bdot)

!      Written 15 January 1999. jlm
!      modified 16 January 2014 (f90). jlm
!      copyright INESC-MN/Jose Luis Martins


!      Version 4.94

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       real(REAL64), intent(in)           ::  adot(3,3)                  !< metric in direct space in atomic units (Bohr radius)

!      output

       real(REAL64), intent(out)          ::  vcell                      !< cell volume
       real(REAL64), intent(out)          ::  bdot(3,3)                  !< metric in reciprocal space (contravariant components)

!      counters
       integer       ::  i, j

!      parameters
       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter :: EPS = 0.0000000001_REAL64


!      compute the minors of adot

       bdot(1,1) = adot(2,2)*adot(3,3) - adot(3,2)*adot(2,3)
       bdot(2,1) = adot(3,2)*adot(1,3) - adot(1,2)*adot(3,3)
       bdot(3,1) = adot(1,2)*adot(2,3) - adot(2,2)*adot(1,3)
       bdot(1,2) = adot(2,3)*adot(3,1) - adot(3,3)*adot(2,1)
       bdot(2,2) = adot(3,3)*adot(1,1) - adot(1,3)*adot(3,1)
       bdot(3,2) = adot(1,3)*adot(2,1) - adot(2,3)*adot(1,1)
       bdot(1,3) = adot(2,1)*adot(3,2) - adot(3,1)*adot(2,2)
       bdot(2,3) = adot(3,1)*adot(1,2) - adot(1,1)*adot(3,2)
       bdot(3,3) = adot(1,1)*adot(2,2) - adot(2,1)*adot(1,2)

!      cell volume (squared)

       vcell = bdot(1,1)*adot(1,1) + bdot(2,1)*adot(2,1) +               &
     &         bdot(3,1)*adot(3,1)
       if(vcell < EPS) then
         write(6,'("    STOPPED in adot_to_bdot.    Cell volume ",       &
     &      "squared= ",e12.4)') vcell

         stop

       endif

       do j=1,3
       do i=1,3
         bdot(i,j) = 2*PI*2*PI*bdot(i,j) / vcell
       enddo
       enddo

       vcell = sqrt(vcell)

       return
       end subroutine adot_to_bdot
