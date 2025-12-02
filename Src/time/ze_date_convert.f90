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

!>  Converts the date between two formats (with and without Y2K bug:-)
!>
!>  \author       José Luís Martins
!>  \version      5.12
!>  \date         2 December 2025.
!>  \copyright    GNU Public License v2

subroutine zedate_convert(dmydate, ymddate, ly2k)

! Written 2 December2025 or is it 2025.12.02 based on ze_date. JLM

  implicit none

! input and output

  logical, intent(in)                ::  ly2k                              !<  if true outputs date with Y2K bug:-)

  character(len=9), intent(inout)    ::  dmydate                           !<  date in the day-month-year format
  character(len=10), intent(inout)   ::  ymddate                           !<  date in the year.month.day format


! local variables

  integer            ::  nm
  character(len=3)   ::  month(13)
  integer            ::  ny
  character(len=2)   ::  cent

! counter

  integer    ::  n

  data month/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Unk'/



  if(ly2k) then

    nm = 13

    if(ymddate(6:7) == '01') then
      nm = 1
    elseif(ymddate(6:7) == '02') then
      nm = 2
    elseif(ymddate(6:7) == '03') then
      nm = 3
    elseif(ymddate(6:7) == '04') then
      nm = 4
    elseif(ymddate(6:7) == '05') then
      nm = 5
    elseif(ymddate(6:7) == '06') then
      nm = 6
    elseif(ymddate(6:7) == '07') then
      nm = 7
    elseif(ymddate(6:7) == '08') then
      nm = 8
    elseif(ymddate(6:7) == '09') then
      nm = 9
    elseif(ymddate(6:7) == '10') then
      nm = 10
    elseif(ymddate(6:7) == '11') then
      nm = 11
    elseif(ymddate(6:7) == '12') then
      nm = 12
    endif

    write(dmydate,'(a2,"-",a3,"-",a2)') ymddate(9:10),month(nm),ymddate(3:4)

  else

    read(dmydate(8:9),'(i2)') ny

    if(ny > 50 ) then
      cent = '19'
    else
      cent = '20'
    endif

    write(ymddate, '(a2,a2,".--.",a2)') cent, dmydate(8:9), dmydate(1:2)

    do n = 1,12

      if(dmydate(4:6) == month(n)) then

        if(n < 10) then
          write(ymddate, '(a2,a2,".0",i1,".",a2)') cent, dmydate(8:9), n, dmydate(1:2)
        else
          write(ymddate, '(a2,a2,".",i2,".",a2)') cent, dmydate(8:9), n, dmydate(1:2)
        endif

        exit

      endif

    enddo

  endif

  return

end subroutine zedate_convert
