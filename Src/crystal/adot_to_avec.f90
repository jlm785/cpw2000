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

!>     given a metric adot, returns an upper-triangular matrix
!>     avec={a1,a2,a3}  such that adot=transpose(avec)*avec

       subroutine adot_to_avec(adot,avec,bvec)

!      Written by Ivo Souza
!      Modified 15 January 1999. jlm
!      Converted to f90. JLM
!      copyright inesc-mn/Jose Luis Martins/Ivo Souza

!      Version 4.94 of cpw

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       real(REAL64), intent(in)           ::  adot(3,3)                  !< metric in direct space in atomic units (Bohr radius)

!      output

       real(REAL64), intent(out)          ::  avec(3,3)                  !< primitive lattice vectors that generate adot in canonical orientation
       real(REAL64), intent(out)          ::  bvec(3,3)                  !< reciprocal lattice vectors

!      parameters

       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter :: EPS = 0.0000000001_REAL64
       real(REAL64), parameter :: ZERO = 0.0_REAL64


!      compute a set of cell vectors compatible with the given metric

       avec(1,1) = adot(1,1)
       if(avec(1,1) < eps) then
         write(6,*)
         write(6,*) '   *** STOPPED in adot_to_avec'
         write(6,'("   avec(1,1) = ",e14.6)') avec(1,1)

         stop

       endif

       avec(1,1) = sqrt(avec(1,1))
       avec(2,1) = ZERO
       avec(3,1) = ZERO

       avec(1,2) = adot(1,2)/avec(1,1)
       avec(2,2) = adot(2,2)-avec(1,2)*avec(1,2)
       if(avec(2,2) < eps) then
         write(6,*)
         write(6,*) '   *** STOPPED in adot_to_avec'
         write(6,'("   avec(2,2) = ",e14.6)') avec(2,2)

         stop

       endif
       avec(2,2) = sqrt(avec(2,2))
       avec(3,2) = ZERO

       avec(1,3) = adot(1,3)/avec(1,1)
       avec(2,3) = (adot(2,3)-avec(1,2)*avec(1,3))/avec(2,2)
       avec(3,3) = adot(3,3) - avec(2,3)*avec(2,3) -                     &
     &                         avec(1,3)*avec(1,3)
       if(avec(3,3) < eps) then
         write(6,*)
         write(6,*) '   *** STOPPED in adot_to_avec'
         write(6,'("   avec(3,3) = ",e14.6)') avec(3,3)

         stop

       endif  
       avec(3,3) = sqrt(avec(3,3))

       bvec(1,1) = 2*pi / avec(1,1)      
       bvec(2,1) =-2*pi*avec(1,2) / (avec(1,1)*avec(2,2))
       bvec(3,1) = 2*pi*avec(1,2)*avec(2,3) /                            &
     &                (avec(1,1)*avec(2,2)*avec(3,3)) -                  &
     &             2*pi*avec(1,3) / (avec(1,1)*avec(3,3))
       bvec(1,2) = ZERO      
       bvec(2,2) = 2*pi / avec(2,2)
       bvec(3,2) =-2*pi*avec(2,3) / (avec(2,2)*avec(3,3))
       bvec(1,3) = ZERO      
       bvec(2,3) = ZERO
       bvec(3,3) = 2*pi / avec(3,3)

       return
       end subroutine adot_to_avec
