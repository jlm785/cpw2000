!------------------------------------------------------------!
!                                                            !
! Copyright (C) 2001 PWSCF group                             !
! This file is distributed under the terms of the            !
! GNU General Public License. See the file `License'         !
! in the root directory of the present distribution,         !
! or http://www.gnu.org/copyleft/gpl.txt .                   !
!                                                            !
! adapted 22 November 2022.  Jose Luis Martins               !
!                                                            !
!------------------------------------------------------------!

!> Transforms a COMPLEX tensor (like the dynamical matrix) from
!> crystal to cartesian axis (\(\text{iflg}\geq 1\)) or viceversa
!> (\(\text{iflg} \leq -1\)).
!>
!>  \author       Marsamos, Fabrizio2 Quantum Espresso, Adapted by Jose Luis Martins
!>  \version      5.06
!>  \date         21 September 2011, 22 November 2022.
!>  \copyright    GNU Public License v2

! adapted by Jose Luis Martins, INESC MN, 22 November 2022.

subroutine QE_trntnsc (phi, avec, bvec, iflg)
!
!  USE kinds, only : DP
!
  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

  integer, intent(in)                ::  iflg                            !<  if > 0: lattice to cartesian; else:  cartesian to lattice
  real(REAL64), intent(in)           ::  avec(3,3)                       !<  primitive lattice vectors in canonical orientation
  real(REAL64), intent(in)           ::  bvec(3,3)                       !<  reciprocal lattice vectors

  complex(REAL64), intent(inout)     ::  phi(3,3)                        !<  the matrix to transform

! local arrays

  complex(REAL64) :: wrk (3, 3)

! parameters

  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  ZERO = 0.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer :: i, j, k, l



  if (iflg > 0) then

!   forward transformation (crystal to cartesian axis)

!    call zcopy (9, phi, 1, wrk, 1)
     wrk(:,:) = phi(:,:)

     do i = 1, 3
        do j = 1, 3
           phi (i, j) = C_ZERO

           do k = 1, 3
              do l = 1, 3
                 phi (i, j) = phi (i, j) + bvec (i, k) * wrk (k, l) * bvec (j, l)
              enddo
           enddo

           phi (i, j) = phi (i, j) / (4*PI*PI)

        enddo
     enddo

  else

!   backward transformation (cartesian to crystal axis)

     do i = 1, 3
        do j = 1, 3
           wrk (i, j) = C_ZERO

           do k = 1, 3
              do l = 1, 3
                 wrk (i, j) = wrk (i, j) + avec (k, i) * phi (k, l) * avec (l, j)
              enddo
           enddo

        enddo
     enddo

!    call zcopy (9, wrk, 1, phi, 1)
     phi(:,:) = wrk(:,:)

  endif

  return

end subroutine QE_trntnsc

