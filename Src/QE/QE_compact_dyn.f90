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

!> writes the dynamical matrix from a 3,3,nat,nat array to a 3*nat,3*nat array
!>
!>  \author       Dal Corso Quantum Espresso, Adapted by Jose Luis Martins
!>  \version      5.06
!>  \date         19 January 2013, 4 December 2022.
!>  \copyright    GNU Public License v2

! adapted by Jose Luis Martins, INESC MN, 22 November 2022

subroutine QE_compact_dyn(nat, dyn, phi)

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

  integer, intent(in)                ::  nat
  complex(REAL64), intent(in)        ::  phi(3,3,nat,nat)
  complex(REAL64), intent(out)       ::  dyn(3*nat, 3*nat)

  integer :: na, nb, icart, jcart, imode, jmode

  do na = 1, nat
     do icart = 1, 3
        imode = 3 * ( na - 1 ) + icart
        do nb = 1, nat
           do jcart = 1, 3
              jmode = 3 * ( nb - 1 ) + jcart
              dyn (imode, jmode) = phi (icart, jcart, na, nb)
           end do
        end do
     end do
  end do

  return

end subroutine QE_compact_dyn

