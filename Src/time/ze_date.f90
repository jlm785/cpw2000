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

!>     Gets the date (day-month-year)
!>     linux gfortran/ifort version
!>     let's code the y2k bug

       subroutine zedate(bdate)

!      Written 21 october 2003
!      Modified 11 September 2015. f90. JLM
!      Modified, documentation, December 2019. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94

       implicit none

!      output

       character(len=9), intent(out)      ::  bdate                      !<  date the subroutine was called

!      local variables

       integer            ::  nm
       character(len=3)   ::  month(13)
       character(len=8)   ::  date
       character(len=10)  ::  time

       data month/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep', &
     &            'Oct','Nov','Dec','Unk'/

       call date_and_time(date, time )

       nm = 13

       if(date(5:6) == '01') then
         nm = 1
       elseif(date(5:6) == '02') then
         nm = 2
       elseif(date(5:6) == '03') then
         nm = 3
       elseif(date(5:6) == '04') then
         nm = 4
       elseif(date(5:6) == '05') then
         nm = 5
       elseif(date(5:6) == '06') then
         nm = 6
       elseif(date(5:6) == '07') then
         nm = 7
       elseif(date(5:6) == '08') then
         nm = 8
       elseif(date(5:6) == '09') then
         nm = 9
       elseif(date(5:6) == '10') then
         nm = 10
       elseif(date(5:6) == '11') then
         nm = 11
       elseif(date(5:6) == '12') then
         nm = 12
       endif

       write(bdate,'(a2,"-",a3,"-",a2)') date(7:8),month(nm),date(3:4)
       
       return
       end subroutine zedate
