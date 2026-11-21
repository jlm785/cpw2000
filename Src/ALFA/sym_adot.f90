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

!>  Symmetrizes the metric in direct space with the
!>  operations of the space group.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         19 November 2026.
!>  \copyright    GNU Public License v2

subroutine sym_adot(adot, ntrans, mtrx)

! Written 19 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

! input and output

  real(REAL64), intent(inout)        ::  adot(3,3)                       !<  metric in direct space

! local variables

  real(REAL64)       ::  adsym(3,3)

! parameters

  real(REAL64), parameter :: ZERO = 0.0_REAL64

! counters

  integer       ::  i, j, k, l, n


! mtrx acts on reciprocal lattice coordinates.  The same operation in direct
! lattice coordinates is  d = (mtrx^-1)^T,  as G.r is invariant.  In lattice
! coordinates mtrx is not orthogonal,  mtrx^-1 /= mtrx^T  in general.
! The metric is invariant if  d^T adot d = mtrx^-1 adot (mtrx^-1)^T = adot.
! The inverses are also in the group, so the sum can be done with mtrx.

  do j = 1,3
  do i = 1,3
    adsym(i,j) = ZERO
  enddo
  enddo

  do n = 1,ntrans
    do j = 1,3
    do i = 1,3
      do l = 1,3
      do k = 1,3
        adsym(i,j) = adsym(i,j) + mtrx(i,k,n)*adot(k,l)*mtrx(j,l,n)
      enddo
      enddo
    enddo
    enddo
  enddo

  do j = 1,3
  do i = 1,3
    adot(i,j) = adsym(i,j) / ntrans
  enddo
  enddo

  return

end subroutine sym_adot
