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

       subroutine zesec(tback)

!      Gets cpu time in seconds
!      Dummy version

!      Written 21 october 2003
!      Modified 11 September 2015. f90. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.60

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      output

       real(REAL64), intent(out)          ::  tback                    !  cpu time in seconds

       tback = 0.0

       return
       end
