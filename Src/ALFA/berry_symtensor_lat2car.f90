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

!>  Converts a symmetric tensor from lattice to cartesian
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         24 January 2023.
!>  \copyright    GNU Public License v2

subroutine berry_symtensor_lat2car(adot, stensor, stensor_car, pv, paxis)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space
  real(REAL64), intent(in)           ::  stensor(3,3)                    !<  Symmetric tensor  (lattice coordinates)

! output

  real(REAL64), intent(out)          ::  stensor_car(3,3)                !<  Symmetric tensor (cartesian coordinates)
  real(REAL64), intent(out)          ::  pv(3)                           !<  principal values
  real(REAL64), intent(out)          ::  paxis(3,3)                      !<  principal axis


! local variable
  real(REAL64)      ::  avec(3,3)           !  primitive lattice vectors
  real(REAL64)      ::  bvec(3,3)           !  reciprocal primitive lattice vectors

  integer           ::  info                !  error from lapack

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64
  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64

! counters

  integer    ::  j, k, m, i


  call adot_to_avec_sym(adot, avec, bvec)

  do j = 1,3
  do k = 1,3
    stensor_car(j,k) = ZERO
    do m = 1,3
    do i = 1,3
      stensor_car(j,k) = stensor_car(j,k) + stensor(m,i)*avec(j,m)*avec(k,i)
    enddo
    enddo
    stensor_car(j,k) = stensor_car(j,k) / (4*PI*PI)
  enddo
  enddo

  call diag_r8(3, stensor_car, pv, paxis, 3, info)

  return

end subroutine berry_symtensor_lat2car
