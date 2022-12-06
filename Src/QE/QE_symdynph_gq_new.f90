!------------------------------------------------------------!
!                                                            !
! Copyright (C) 2001-2012 Quantum ESPRESSO group             !
! This file is distributed under the terms of the            !
! GNU General Public License. See the file `License'         !
! in the root directory of the present distribution,         !
! or http://www.gnu.org/copyleft/gpl.txt .                   !
!                                                            !
! adapted 2 December 2022.  Jose Luis Martins                !
!                                                            !
!------------------------------------------------------------!

!>  This routine receives as input an unsymmetrized dynamical
!>  matrix expressed on the crystal axes and imposes the symmetry
!>  of the small group of q. Furthermore it imposes also the symmetry
!>  q -> -q+G if present.
!>  February 2020: Update (A. Urru) to include the symmetry operations
!>  that require the time reversal operator (meaning that TS is a
!>  symmetry of the crystal). For more information please see:
!>  Phys. Rev. B 100, 045115 (2019).
!>
!>  \author       marsamos,...,fabrizio2, Quantum Espresso, Adapted by Jose Luis Martins
!>  \version      5.06
!>  \date         21 September 2011, 9 September 2021.
!>  \copyright    GNU Public License v2





subroutine QE_symdynph_gq_new( xq, phi, s, invs, rtau, irt, nsymq,       &
                            nat, irotmq, minus_q, t_rev)

!   USE kinds, only : DP
!   USE constants, ONLY: tpi
!   USE symm_base, ONLY : t_rev

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  xq(3)                           !< the q point, lattice coordinates

  integer, intent(in)                ::  nat                             !<  the number of atoms

  integer, intent(in)                ::  nsymq                           !<  the order of the small group
  integer, intent(in)                ::  s(3,3,48)                       !<  the symmetry matrices
  integer, intent(in)                ::  invs(48)                        !<  the inverse of each matrix

  integer, intent(in)                ::  irt(48,nat)                     !<  the rotated of each atomic vector
  integer, intent(in)                ::  irotmq                          !<  the rotation sending q ->-q+G

  real(REAL64), intent(in)           ::  rtau(3,48,nat)                  !<  rotated atom minus original, lattice coordinates.

  logical, intent(in)                ::  minus_q                         !<  true if a symmetry q->-q+G is present
  integer, intent(in)                ::  t_rev(48)                       !<  1 if symmetry includes time reversal

! input and output

  complex(REAL64), intent(inout)     ::  phi(3,3,nat,nat)                !<  the matrix to symmetrize

! local variables

  integer            ::  sna, snb, irot, iflb (nat, nat)                 ! indices, work space

  real(REAL64)       ::  arg                                             ! the argument of the phase

  complex(REAL64)    ::  phip (3, 3, nat, nat)                           ! work space
  complex(REAL64)    ::  work (3, 3), fase, faseq (48)                   ! work space, phase factors

! parameters

  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  ZERO = 0.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer :: isymq, na, nb, ipol, jpol, lpol, kpol

!    We start by imposing hermiticity

  do na = 1, nat
     do nb = 1, nat
        do ipol = 1, 3
           do jpol = 1, 3
              phi (ipol, jpol, na, nb) = ( phi (ipol, jpol, na, nb)      &
                   + CONJG(phi (jpol, ipol, nb, na) ) ) / 2
              phi (jpol, ipol, nb, na) = CONJG(phi (ipol, jpol, na, nb) )
           enddo
        enddo
     enddo
  enddo

! If no other symmetry is present we quit here

  if ( (nsymq == 1) .and. (.not.minus_q) ) RETURN

! Then we impose the symmetry q -> -q+G if present

  if (minus_q) then
     do na = 1, nat
        do nb = 1, nat
           do ipol = 1, 3
              do jpol = 1, 3
                 work(:,:) = C_ZERO
                 sna = irt (irotmq, na)
                 snb = irt (irotmq, nb)
                 arg = ZERO
                 do kpol = 1, 3
                    arg = arg + (xq (kpol) * (rtau (kpol, irotmq, na) -  &
                                              rtau (kpol, irotmq, nb) ) )
                 enddo
                 arg = arg * 2*PI
                 fase = CMPLX(cos (arg), sin (arg) , REAL64)
                 do kpol = 1, 3
                    do lpol = 1, 3
                       work (ipol, jpol) = work (ipol, jpol) +                &
                            s (ipol, kpol, irotmq) * s (jpol, lpol, irotmq)   &
                            * phi (kpol, lpol, sna, snb) * fase
                    enddo
                 enddo
                 phip (ipol, jpol, na, nb) = (phi (ipol, jpol, na, nb) +      &
                      CONJG( work (ipol, jpol) ) ) / 2
              enddo
           enddo
        enddo
     enddo
     phi = phip
  endif


! Here we symmetrize with respect to the small group of q

  if (nsymq == 1) RETURN

  iflb (:, :) = 0
  do na = 1, nat
     do nb = 1, nat
        if (iflb (na, nb) == 0) then
           work(:,:) = C_ZERO
           do isymq = 1, nsymq
              irot = isymq
              sna = irt (irot, na)
              snb = irt (irot, nb)
              arg = ZERO
              do ipol = 1, 3
                 arg = arg + (xq (ipol) * (rtau (ipol, irot, na) -       &
                                           rtau (ipol, irot, nb) ) )
              enddo
              arg = arg * 2+PI
              faseq (isymq) = CMPLX(cos (arg), sin (arg) , REAL64)
              do ipol = 1, 3
                 do jpol = 1, 3
                    do kpol = 1, 3
                       do lpol = 1, 3
                          IF (t_rev(isymq) == 1) THEN
                             work (ipol, jpol) = work (ipol, jpol) +            &
                                  s (ipol, kpol, irot) * s (jpol, lpol, irot)   &
                           * CONJG(phi (kpol, lpol, sna, snb) * faseq (isymq))
                          ELSE
                             work (ipol, jpol) = work (ipol, jpol) +            &
                                  s (ipol, kpol, irot) * s (jpol, lpol, irot)   &
                                 * phi (kpol, lpol, sna, snb) * faseq (isymq)
                          ENDIF
                       enddo
                    enddo
                 enddo
              enddo
           enddo
           do isymq = 1, nsymq
              irot = isymq
              sna = irt (irot, na)
              snb = irt (irot, nb)
              do ipol = 1, 3
                 do jpol = 1, 3
                    phi (ipol, jpol, sna, snb) = C_ZERO
                    do kpol = 1, 3
                       do lpol = 1, 3
                          IF (t_rev(isymq)==1) THEN
                             phi(ipol,jpol,sna,snb) = phi(ipol,jpol,sna,snb)     &
                             + s(ipol,kpol,invs(irot))*s(jpol,lpol,invs(irot))   &
                               * CONJG(work (kpol, lpol)*faseq (isymq))
                          ELSE
                             phi(ipol,jpol,sna,snb) = phi(ipol,jpol,sna,snb)     &
                             + s(ipol,kpol,invs(irot))*s(jpol,lpol,invs(irot))   &
                               * work (kpol, lpol) * CONJG(faseq (isymq) )
                          ENDIF
                       enddo
                    enddo
                 enddo
              enddo
              iflb (sna, snb) = 1
           enddo
        endif
     enddo
  enddo

  phi (:, :, :, :) = phi (:, :, :, :) / nsymq

  return

end subroutine QE_symdynph_gq_new

