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

!>     calculates the Fermi-Dirac function
!>     fd = 1/(exp(e)+1) and its derivative
!>     15 digit accuracy tested against Mathematica

       subroutine fermi_dirac(en,fd,dfdden)

!      written 17 June 2014.
!      Modified documentation august 2019.  JLM
!      copyright Jose Luis Martins / INESC-MN

!      version 4.94

       implicit none
       integer, parameter  :: REAL64 = selected_real_kind(12)

!      input

       real(REAL64), intent(in)         ::  en                           !<  energy minus chemical potential

!      output

       real(REAL64), intent(out)        ::  fd                           !<  Fermi-Dirac function 1/(exp(e)+1)
       real(REAL64), intent(out)        ::  dfdden                       !<  d fd / d en

!      constants

       real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

       if(abs(en) < 25) then
         fd = UM / (exp(en) + UM)
         dfdden = -fd*fd * exp(en)
       else
         if(en > ZERO) then
           fd = exp(-en)*(UM-exp(-en))
           dfdden = -exp(-en)*(UM-2*exp(-en))
         else
           fd = UM - exp(en)
           dfdden = - exp(en)*(UM-2*exp(en))
         endif
       endif

       return

       end
