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

!>  Weights for a quadratic interpolation on a 3D grid using 10 points.
!>  Notice that it is asymetric. Priviliged direction [111]

subroutine quad_3D10pt(d,fg)

! Written 29 November 2013, based on the old plotting code. JLM 
! Modified documentation August 2019.  JLM
! Modified to output weights, 17 October 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.98

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)        ::  d(3)                               !<  point to be interpolated in grid. Should be 0 < d(i) < 1.

! output

  real(REAL64), intent(out)       ::  fg(-1:1,-1:1,-1:1)                 !<  weights on a 3x3x3 grid

! constants

  real(REAL64), parameter ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64


  fg(:,:,:) = ZERO

  fg( 0, 0, 0) = UM + d(1)*d(2) + d(2)*d(3) + d(3)*d(1) - d(1)*d(1) - d(2)*d(2) - d(3)*d(3)

  fg(-1, 0, 0) = d(1) * ( d(1) - UM ) / 2
  fg( 0,-1, 0) = d(2) * ( d(2) - UM ) / 2
  fg( 0, 0,-1) = d(3) * ( d(3) - UM ) / 2

  fg( 1, 0, 0) = d(1) * ( UM + d(1) - 2*d(2) -2*d(3) ) / 2
  fg( 0, 1, 0) = d(2) * ( UM + d(2) - 2*d(3) -2*d(1) ) / 2
  fg( 0, 0, 1) = d(3) * ( UM + d(3) - 2*d(1) -2*d(2) ) / 2

  fg( 0, 1, 1) = d(2) * d(3)
  fg( 1, 0, 1) = d(3) * d(1)
  fg( 1, 1, 0) = d(1) * d(2)

  return
end subroutine quad_3D10pt
