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

!>  a non-advancing status counter.
!>
!>  \author       Carlos Loia Reis
!>  \version      5.13
!>  \date         Unknown, 4 October 2026.
!>  \copyright    GNU Public License v2


subroutine progress(j,n)

! Implemented by Carlos Loia Reis at an unknown date.
! Documentation 5 october 2026. JLM

  implicit none

! input

  integer, intent(in)                   ::  j                            !<  current step
  integer, intent(in)                   ::  n                            !<  number of steps

  !  note that achar(13) brings the cursor to the begining of line

  write(6,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &
  & " Percent Complete: ", (real(j)/real(n))*100.0, "%"

  return

end subroutine progress
