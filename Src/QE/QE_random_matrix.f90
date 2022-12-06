!--------------------------------------------------------------------------------------!
!                                                                                      !
! Copyright (C) 2001 PWSCF group                                                       !
! This file is distributed under the terms of the                                      !
! GNU General Public License. See the file `License'                                   !
! in the root directory of the present distribution,                                   !
! or http://www.gnu.org/copyleft/gpl.txt .                                             !
!                                                                                      !
! The random matrix (see later) can be populated either with uniformly                 !
! distributed random numbers or with normal-distributed random numbers.                !
! The former has been the default until QE 6.1, however it sometimes                   !
! produces accidentally degenerate eigenvalue, especially when dealing                 !
! with large number of atoms.                                                          !
! A matrix of normal-distributed numbers should have (on average) more                 !
! evenly spaced eigenvalues, reducing the chance of collision.                         !
!                                                                                      !
! See <http://web.math.princeton.edu/mathlab/projects/ranmatrices/yl/randmtx.PDF>      !
! (If I understand it correctly)                                                       !
! LP 2017                                                                              !
! Uncomment the following line in case of trouble with set_irr_sym_new                 !
!!#define __UNIFORM_DISTRIB                                                            !
!                                                                                      !
! adapted 21 December 2022.  Jose Luis Martins                                         !
!                                                                                      !
!--------------------------------------------------------------------------------------!

!> Create a Random hermitian matrix with non zero elements similar to
!> the dynamical matrix of the system.
!>
!>  \author       Marsamos, Fabrizio2 Quantum Espresso, Adapted by Jose Luis Martins
!>  \version      5.06
!>  \date         21 September 2011, 22 November 2022.
!>  \copyright    GNU Public License v2

subroutine QE_random_matrix_new (nat, irt, nsymq, lminus_q, irotmq,      &
     wdyn, lgamma, idum)

! #if defined (__UNIFORM_DISTRIB)
! #define __RANDOM_DBLE  CMPLX(2.0_DP*randy () - 1.0_DP, 0.d0,kind=DP)
! #define __RANDOM_CMPLX CMPLX(2.0_DP * randy () - 1.0_DP, 2.0_DP * randy () - 1.0_DP,kind=DP)
! #else
! #define __RANDOM_DBLE  gauss_dist_scal(0._dp, 1._dp)
! #define __RANDOM_CMPLX gauss_dist_cmplx(0._dp, 1._dp)
!   USE random_numbers, ONLY : gauss_dist_scal, gauss_dist_cmplx
! #endif

!   USE kinds, only : DP
!   USE random_numbers, ONLY : randy

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  nat                             !<  Total number of atoms
  integer, intent(in)                ::  irt(48,nat)                     !<  index of the rotated atom
  integer, intent(in)                ::  nsymq                           !<  the order of the small group
  integer, intent(in)                ::  irotmq                          !<  the rotation sending q -> -q
  logical, intent(in)                ::  lgamma                          !<  if TRUE  q=0
  logical, intent(in)                ::  lminus_q                        !<  if TRUE there is an inversion symmetry

  integer, intent(inout)             ::  idum                            !<  random number seed

! output

  complex(REAL64), intent(out)       :: wdyn (3, 3, nat, nat)            !<  random matrix

! local variables

  integer         ::  ira                       !  rotated atom
  integer         ::  iramq                     !  rotated atom with the q->-q+G symmetry

  real(REAL64)    ::  rr, rc

  real(REAL64),external      :: ran2

! parameters

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer ::  na, nb, ipol, jpol, isymq




  wdyn (:, :, :, :) = C_ZERO

  do na = 1, nat
     do ipol = 1, 3

        rr = 2*ran2(idum) - UM
        wdyn (ipol, ipol, na, na) = 2*cmplx(rr,ZERO,REAL64)

        if(ipol < 3) then
           do jpol = ipol + 1, 3
              rr = 2*ran2(idum) - UM
              if (lgamma) then
                 rc = ZERO
              else
                 rc = 2*ran2(idum) - UM
              endif
              wdyn (ipol, jpol, na, na) = cmplx(rr, rc,REAL64)
              wdyn (jpol, ipol, na, na) = cmplx(rr,-rc,REAL64)
           enddo
        endif

        if(na < nat) then
           do nb = na + 1, nat

              do isymq = 1, nsymq
                 ira = irt (isymq, na)
                 if (lminus_q) then
                    iramq = irt (irotmq, na)
                 else
                    iramq = 0
                 endif

                 if ( (nb == ira) .or. (nb == iramq) ) then
                    do jpol = 1, 3
                       rr = 2*ran2(idum) - UM
                       if (lgamma) then
                          rc = ZERO
                       else
                          rc = 2*ran2(idum) - UM
                       endif
                       wdyn(ipol, jpol, na, nb) = cmplx(rr, rc,REAL64)
                       wdyn(jpol, ipol, nb, na) = cmplx(rr,-rc,REAL64)
                    enddo

                    exit

                 endif

              enddo


           enddo
        endif

     enddo
  enddo

  return

end subroutine QE_random_matrix_new
