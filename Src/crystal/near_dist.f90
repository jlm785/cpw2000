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

!>  Closest distance between atoms taking into account periodicity.
!>  Calls near_dist_image (exact for any cell, the distance is not squared).
!>
!>  \author       José Luís Martins
!>  \version      5.13
!>  \date         30 November 2016, 10 October 2026.
!>  \copyright    GNU Public License v2

subroutine near_dist(distmin, adot, r0, r1)

! Written 30 November 2016.  JLM
! Modified, documentation, August 2019. JLM
! Documentation, missing or incomplete argument description. 28 September 2026. JLM+claude
! Indentation. 28 September 2026. JLM+claude
! Calls near_dist_image: exact for any cell, returns the distance (was its square). 10 October 2026. JLM+claude


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  real(REAL64), intent(in)           ::  r0(3)                           !<  position of the first atom (lattice coordinates)
  real(REAL64), intent(in)           ::  r1(3)                           !<  position of the second atom (lattice coordinates)

! output

  real(REAL64), intent(out)          ::  distmin                         !<  closest distance between images (not squared)

! local variables

  real(REAL64)    ::  rmin(3)                                            !  nearest image of r1 - r0 (not used)


  call near_dist_image(distmin, rmin, adot, r0, r1)

  return

end subroutine near_dist

