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

!>     gets the elapsed time since beginning of month
!>     precision of miliseconds
!>     linux gfortran/ifort version

       subroutine zeelap(el_time)

!      Written 2 may 2006
!      Modified 11 September 2015. f90. JLM
!      Modified, documentation, beginning of month, December 2019. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      output

       real(REAL64), intent(out)          ::  el_time                    !<  elapsed time since beginning of month

!      local variables

       character(len=8)   ::  date
       character(len=10)  ::  time
       character(len=5)   ::  zone
       integer            ::  ival(8)

!      constants

       real(REAL64), parameter  ::  UM = 1.0_REAL64

       call date_and_time(date,time,zone,ival)

       el_time = ((((ival(3)-1)*24 + ival(5))*60 + ival(6))*60 +         &
     &               ival(7))*UM + (ival(8)*UM)/(1000*UM)

       return
       end subroutine zeelap
