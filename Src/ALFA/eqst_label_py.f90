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

!>  Converts a label for gnuplot into a label for a python string.
!>  The gnuplot escape of the underscore, \_, is replaced by _,
!>  and double quotes are replaced by single quotes.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         4 October 2026.
!>  \copyright    GNU Public License v2

subroutine eqst_label_py(label, pylabel)

! Written 4 October 2026. JLM+claude

  implicit none

! input

  character(len=*), intent(in)          ::  label                        !<  label for gnuplot

! output

  character(len=*), intent(out)         ::  pylabel                      !<  label for python

! local variables

  integer                   ::  i, k, n


  pylabel = ' '
  n = len_trim(label)

  i = 1
  k = 0
  do while(i <= n .and. k < len(pylabel))
    k = k + 1
    if(i < n .and. label(i:min(i+1,n)) == '\_') then
      pylabel(k:k) = '_'
      i = i + 2
    elseif(label(i:i) == '"') then
      pylabel(k:k) = "'"
      i = i + 1
    else
      pylabel(k:k) = label(i:i)
      i = i + 1
    endif
  enddo

  return

end subroutine eqst_label_py
