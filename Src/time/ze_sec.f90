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

!>     Gets cpu time in seconds (sum over cores)
!>     linux gfortran/ifort version

       subroutine zesec(tback)

!      Gets cpu time in seconds
!      linux gfortran/ifort version

!      Written 21 october 2003
!      Modified 11 September 2015. f90. JLM
!      Modified, documentation, December 2019. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)
       integer, parameter          :: REAL32 = selected_real_kind(6)

!      output

       real(REAL64), intent(out)          ::  tback                    !<  cpu time in seconds

!      local variables

       real(REAL32)    :: tin

       call cpu_time(tin)
       tback = tin

       return
       end subroutine zesec
