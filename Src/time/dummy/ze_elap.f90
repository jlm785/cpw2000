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

!>  Gets the elapsed time since midnight.
!>  Dummy version.
!>
!>  \author       José Luís Martins
!>  \version      5.13
!>  \date         2 May 2006, 11 September 2015.
!>  \copyright    GNU Public License v2

subroutine zeelap(el_time)

! Written 2 may 2006
! Modified 11 September 2015. f90. JLM
! Indentation. 28 September 2026. JLM+claude


  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! output

  real(REAL64), intent(out)          ::  el_time                         !<  dummy elapsed time since midnight

  el_time = 0.0

  return

end subroutine zeelap
