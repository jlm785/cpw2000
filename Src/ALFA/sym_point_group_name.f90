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

!>  Gives the rotation-inversion matrices identifies the point group.
!>  Interface to the Quantum-espresso subroutine
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         23 November 2022.
!>  \copyright    GNU Public License v2

subroutine sym_point_group_name(adot, ipr, code_group,                   &
     ntrans, mtrx)

! Written 23 November 2022.


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ipr                             !<  prints the group name

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

! output

  integer, intent(out)                ::  code_group                     !<  code for the point group (see QE_find_group)

! local variables

  real(REAL64)        ::  tnp(3,48)
  real(REAL64)        ::  rot(3,3,48)
  real(REAL64)        ::  tau(3,48)
  character (len=11)  ::  gname

! parameter

  real(REAL64), parameter :: ZERO = 0.0_REAL64

  tnp(:,:) = ZERO

  call sym_cartesian_op(adot, rot, tau, ntrans, mtrx, tnp)

  call QE_find_group( ntrans, rot, gname, code_group )

  if(ipr > 0) then
    write(6,*)
    write(6,*) '   The point group is:'
    write(6,'(3x,a11)') gname
    write(6,*)
  endif

  return

end subroutine sym_point_group_name
