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

!>  Rotates a vibrational mode
!>
!>  \author       Marsamos, Giannozzi, Dal Corso Quantum Espresso, Adapted by Jose Luis Martins
!>  \version      5.06
!>  \date         21 September 2011, 9 September 2021.
!>  \copyright    GNU Public License v2

! adapted by Jose Luis Martins, INESC MN, 5 December 2022.

! instead of rotating all modes by one symmetry, it rotates one mode by all symmetries.

 subroutine QE_rotate_mod(xq, nat, mode, rmode,                          &
        ntrans, sr, irt, rtau, invs)

!   use kinds, only : dp
!   use constants, only: tpi

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  nat                             !<  the number of atoms

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the little group

  real(REAL64), intent(in)           ::  xq(3)                           !< the q vector

  integer, intent(in)                ::  invs(48)                        !< the index of the inverse of the rotation sr
  integer, intent(in)                ::  irt(48,nat)                     !<  the rotated of each atomic vector
  real(REAL64), intent(in)           ::  sr(3,3,48)                      !<  the symmetry matrices in cartesian coordinates
  real(REAL64), intent(in)           ::  rtau(3,48,nat)                  !< the vector r = s tau - tau for all rotations

  complex(REAL64), intent(in)        ::  mode(3*nat)                     !<  the mode to rotate

! output

  complex(REAL64), intent(out)       ::  rmode(3*nat,ntrans)             !<  the rotated mode

! local variables

  complex(REAL64)     ::  phase   ! auxiliary phase
  real(REAL64)        ::  arg     ! an auxiliary argument
  integer             ::  irot

! parameters

  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  ZERO = 0.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters


  integer         ::  n, na, nb, ipol, jpol, mu_i, mu_j


! initialization

  rmode = C_ZERO

  do n = 1,ntrans
     irot = invs(n)

     do na = 1,nat

        nb = irt(irot,na)
        arg = ( xq(1)*rtau(1,irot,na) + xq(2)*rtau(2,irot,na) +  &
                xq(3)*rtau(3,irot,na) ) * (2*PI)
        phase = cmplx(cos(arg), sin(arg), REAL64)

        do ipol = 1,3
           mu_i = 3*(na-1)+ipol
           do jpol = 1,3
              mu_j = 3*(nb-1)+jpol
              rmode(mu_i,n)=rmode(mu_i,n) + sr(ipol,jpol,n)*mode(mu_j)*phase
           enddo
        enddo

     enddo

  enddo

  return

end subroutine QE_rotate_mod

