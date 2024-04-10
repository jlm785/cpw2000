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

!>  Performs an approximate Newton step.
!>  No relaxation on "coarse" grid, jacobi relaxation for fine grid.
!>
!>  Some people call it pre-conditioning...
!>
!>  \author       Jose Luis Martins
!>  \version      5.09
!>  \date         October 4, 1989, January 2020.
!>  \copyright    GNU Public License v2

subroutine diag_rq_jac_c16(phi, lambda, hdiag, mtxd, neig,               &
     mxddim, mxdbnd)


! Written October 4, 1989. JLM
! modified April 1999. JLM (order n**2)
! modified April 22, 2014, f90. JLM
! Modified, documentation, omp, January 2020. JLM
! Modified indentation. 13 December 2023.
! Prefix diag_, 17 March 2024. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  real(REAL64), intent(in)           ::  lambda(mxdbnd)                  !<  expectation values (Hartree)
  real(REAL64), intent(in)           ::  hdiag(mxddim)                   !<  hamiltonian diagonal
  integer, intent(in)                ::  neig                            !<  number of eigenvectors
  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian

! input and output

  complex(REAL64), intent(inout)     ::  phi(mxddim,mxdbnd)              !<  error vectors on input, correction on output

! local variables

  real(REAL64)   ::  x

! constants

  real(REAL64), parameter  :: UM = 1.0_REAL64

! counters

  integer         ::  i,j


  do i = 1,neig
!$omp parallel do default(shared) private(j,x)
    do j = 1,mtxd
      x = hdiag(j) - lambda(i)
      if(x < UM) then
        x = UM
      else
        x = UM/x
      endif
      phi(j,i) = -x * phi(j,i)
    enddo
!$omp end parallel do
  enddo

  return

end subroutine diag_rq_jac_c16
