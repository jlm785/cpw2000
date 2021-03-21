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

       subroutine zeelap(el_time)

!      gets the elapsed time since midnight
!      Dummy version

!      Written 2 may 2006
!      Modified 11 September 2015. f90. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.60

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      output

       real(REAL64), intent(out)          ::  el_time                    !  elapsed time since midnight

       el_time = 0.0

       return
       end
