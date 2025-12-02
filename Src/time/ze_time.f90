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

!>  gets the time of day (hh:mm:ss)
!>  linux gfortran/ifort version
!>
!>  \author       José Luís Martins
!>  \version      5.12
!>  \date         21 october 2003, 2 December 2025.
!>  \copyright    GNU Public License v2

subroutine zetime(btime)

! Written 21 october 2003
! Modified 11 September 2015. f90. JLM
! Modified, documentation, December 2019. JLM
! Modified, indentation, 2 December2025. JLM

  implicit none

! output

  character(len=8), intent(out)      ::  btime                           !<  time of day (hh:mm:ss)

! local variables

  character(len=8)   ::  date
  character(len=10)  ::  time

  call date_and_time(date, time )

  write(btime,'(a2,":",a2,":",a2)') time(1:2),time(3:4),time(5:6)

  return

end subroutine zetime
