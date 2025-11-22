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

!>  Checks if strings have the same info
!>  irrespective of case and leading or trailing blanks.
!>
!>
!>  \author       JLM adapted from the web
!>  \version      5.12
!>  \date         22 November 2025.
!>  \copyright    GNU Public License v2

logical function chrsameinfo(string1, string2)

! Written 22 November 2025. JLM

  implicit none

  character(len=*), intent(in)       ::  string1                         !<  string to be compared
  character(len=*), intent(in)       ::  string2                         !<  string to be compared

! local copies

  character(len=len(string1))        ::  str1
  character(len=len(string2))        ::  str2

  chrsameinfo = .FALSE.

  str1 = adjustl(trim(string1))
  str2 = adjustl(trim(string2))

  call chrcap(str1,0)
  call chrcap(str2,0)

  if(str1 == str2) chrsameinfo = .TRUE.

  return

end function chrsameinfo

