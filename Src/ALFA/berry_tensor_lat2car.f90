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

!>  Converts a tensor of the Berry type from lattice to cartesian
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         24 January 2023, 6 April 2024.
!>  \copyright    GNU Public License v2

subroutine berry_tensor_lat2car(adot, t_lat, t_car, levdegnl, mxddeg)

! adapted from berry_pseudovec_lat2car, 6 April 2024.  JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxddeg                          !<  array dimension for number of levels

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space
  real(REAL64), intent(in)           ::  t_lat(mxddeg,mxddeg,3,3)        !<  Berry tensor in reciprocal lattice coordinates
  integer, intent(in)                ::  levdegnl                        !<  degeneracy of level

! output

  real(REAL64), intent(out)          ::  t_car(mxddeg,mxddeg,3,3)        !<  Berry tensor in cartesian coordinates

! local variables

  real(REAL64)      ::  avec(3,3)           !  primitive lattice vectors
  real(REAL64)      ::  bvec(3,3)           !  reciprocal primitive lattice vectors

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64
  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64

! counters

  integer    ::  i, j, k1, k2, n1, n2


  call adot_to_avec_sym(adot, avec, bvec)

  do n1 = 1,levdegnl
  do n2 = 1,levdegnl

    do k1 = 1,3
    do k2 = 1,3
      t_car(n1,n2,k1,k2) = ZERO
      do j = 1,3
      do i = 1,3
        t_car(n1,n2,k1,k2) = t_car(n1,n2,k1,k2) + avec(k1,i)*t_lat(n1,n2,i,j)*avec(k2,j)
      enddo
      enddo
      t_car(n1,n2,k1,k2) = t_car(n1,n2,k1,k2) / (4*PI*PI)
    enddo
    enddo

  enddo
  enddo

  return

end subroutine berry_tensor_lat2car
