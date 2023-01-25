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

!>  Converts Berry curvature from lattice to cartesian
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         24 January 2023.
!>  \copyright    GNU Public License v2

subroutine berry_curve_2_cart(adot, bcurv, bcurv_car)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space
  real(REAL64), intent(in)           ::  bcurv(3)                        !<  Berry curvature  (lattice coordinates)

! output

  real(REAL64), intent(out)          ::  bcurv_car(3)                    !<  Berry curvature  (cartesian coordinates)

! local variables

  real(REAL64)      ::  tenlat(3,3), tencar(3,3)

  real(REAL64)      ::  avec(3,3)           !  primitive lattice vectors
  real(REAL64)      ::  bvec(3,3)           !  reciprocal primitive lattice vectors

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64
  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64

! counters

  integer    ::  j, k, m, i


  call adot_to_avec_sym(adot, avec, bvec)

  tenlat(1,1) = ZERO
  tenlat(2,2) = ZERO
  tenlat(3,3) = ZERO

  tenlat(2,3) = bcurv(1)
  tenlat(3,1) = bcurv(2)
  tenlat(1,2) = bcurv(3)

  tenlat(3,2) = -tenlat(2,3)
  tenlat(1,3) = -tenlat(3,1)
  tenlat(2,1) = -tenlat(1,2)
  do j = 1,3
  do k = 1,3
    tencar(j,k) = ZERO
    do m = 1,3
    do i = 1,3
      tencar(j,k) = tencar(j,k) + tenlat(m,i)*avec(j,m)*avec(k,i)
    enddo
    enddo
    tencar(j,k) = tencar(j,k) / (4*PI*PI)
  enddo
  enddo

  bcurv_car(1) = tencar(2,3)
  bcurv_car(2) = tencar(3,1)
  bcurv_car(3) = tencar(1,2)

  return

end subroutine berry_curve_2_cart
