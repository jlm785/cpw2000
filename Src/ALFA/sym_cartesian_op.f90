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

!>  Gives the rotation-inversion and fractional translation
!>  in cartesian coordinates (and canonical orientation)
!>  from their lattice coordinate values.
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         23 November 2022.
!>  \copyright    GNU Public License v2

subroutine sym_cartesian_op(adot, rot, tau,                              &
     ntrans, mtrx, tnp)

! Written 23 November 2022.


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(in)           ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

! output

  real(REAL64), intent(out)          ::  rot(3,3,48)                       !<  rotational matrix in cartesian coordinates
  real(REAL64), intent(out)          ::  tau(3,48)                         !<  fractional translation vector in cartesian coordinates

! local variables

  real(REAL64)        ::   avec(3,3), bvec(3,3)
  real(REAL64)        ::   tmp(3,3)

! parameters

  real(REAL64), parameter :: ZERO = 0.0_REAL64
  real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64

! counters

  integer i, j, k, n

! conventional lattice vectors

  call adot_to_avec_sym(adot, avec, bvec)

  do n = 1,ntrans

    do i = 1,3
    do j = 1,3
      tmp(i,j) = ZERO
      do k = 1,3
        tmp(i,j) = tmp(i,j) + bvec(i,k)*mtrx(k,j,n)
      enddo
    enddo
    enddo

    do i = 1,3
    do j = 1,3
      rot(i,j,n) = ZERO
      do k = 1,3
        rot(i,j,n) = rot(i,j,n) + tmp(i,k)*avec(j,k)
      enddo
      rot(i,j,n) = rot(i,j,n) / (2*PI)
    enddo
    enddo

    do i = 1,3
      tau(i,n) = ZERO
      do k = 1,3
        tau(i,n) = tau(i,n) + tnp(k,n)*avec(i,k)
      enddo
      tau(i,n) = tau(i,n) / (2*PI)
    enddo
  enddo

  return

end subroutine sym_cartesian_op
