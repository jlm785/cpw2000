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

!>  Finds improved guess eigen-vector by Jacobian relaxation and
!>  diagonalizes hamiltonian in the initial+improved subspace.
!>  spin version.
!>
!>  \author       Jose Luis Martins
!>  \version      5.09
!>  \date         April 2019. 13 December 2023.
!>  \copyright    GNU Public License v2

subroutine diag_jacobi_ritz_spin_c16(mtxd, neig, eg,                     &
     psi_sp, hpsi_sp, nout,                                              &
     njac, nritz, epsdeg,                                                &
     ng, kgv,                                                            &
     ekpg, isort, vscr_sp, kmscr, nsp,                                   &
     anlsp, xnlkbsp, nanlsp, lnewanl,                                    &
     mxddim, mxdbnd, mxdasp, mxdgve, mxdscr, mxdnsp)

! Adapted from non-spin version, 13 December 2023.


  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !< array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !< array dimension for number of bands
  integer, intent(in)                ::  mxdasp                          !< array dimension of number of projectors
  integer, intent(in)                ::  mxdgve                          !< array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !< array dimension of vscr_sp
  integer, intent(in)                ::  mxdnsp                          !< array dimension for number of spin components (1,2,4)

  integer, intent(in)                ::  mtxd                            !< wavefunction dimension (basis size)
  integer, intent(in)                ::  neig                            !< number of wavefunctions

  integer, intent(in)                ::  njac                            !< number of jacobi iterations
  integer, intent(in)                ::  nritz                           !< number of ritz iterations
  real(REAL64), intent(in)           ::  epsdeg                          !< two eigenvalues are considered quasi degenerate if |e_i - e_j| < eps

  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !< length of k+g-vector of row/column i
  integer, intent(in)                ::  isort(mxddim)                   !< g-vector associated with row/column i of hamiltonian

  real(REAL64), intent(in)           ::  vscr_sp(mxdscr,mxdnsp)          !< screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !< max value of kgv(i,n) used for the potential fft mesh
  integer, intent(in)                ::  nsp                             !< number of spin components ox xc-potential (1,2,4)

  integer, intent(in)                ::  ng                              !< total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !< i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  nanlsp                          !< number of projectors
  complex(REAL64), intent(in)        ::  anlsp(2*mxddim,mxdasp)          !< Kleinman-Bylander projectors
  real(REAL64), intent(in)           ::  xnlkbsp(mxdasp)                 !< Kleinman-Bylander normalization

! output

  complex(REAL64), intent(out)       ::  hpsi_sp(2*mxddim,mxdbnd)        !< H |psi_sp>
  real(REAL64), intent(out)          ::  eg(2*mxdbnd)                    !< <psi_sp| H |psi_sp>
  integer, intent(out)               ::  nout                            !< number of eigenvectors on output  neig <= nout <= mxdbnd, if negative error in the diagonalization.

! input and output

  complex(REAL64), intent(inout)     ::  psi_sp(2*mxddim,mxdbnd)         !< |psi_sp> on input, improved |psi_sp> on output
  logical, intent(inout)             ::  lnewanl                         !< indicates that anlsp has been recalculated (not used in default implementation)

! local variables

  integer           ::  nbas
  integer           ::  mxdorb
  integer           ::  ired

! local allocatable arrays

  complex(REAL64), allocatable       ::  bas(:,:)
  complex(REAL64),allocatable        ::  hbas(:,:)

  integer, allocatable               ::  irow(:)                         !  if irow(i)=0 the vector is linearly dependent on the others and the corresponding xvec is zero

! counters

  integer      ::  j

! work arrays.  Use pointers in case you want to reduce memory requirements...

  allocate(bas(2*mxddim,2*mxdbnd))
  allocate(hbas(2*mxddim,2*mxdbnd))

  do j = 1,nritz

    call zlacpy(' ', 2*mtxd, neig, psi_sp, 2*mxddim, bas, 2*mxddim)

!   Jacobi relaxation

    call diag_jacobi_iter_spin_c16(mtxd, neig, bas, hbas, njac,          &
        bas(1, neig+1), hbas(1, neig+1),                                 &
        ng, kgv,                                                         &
        ekpg, isort, vscr_sp, kmscr, nsp,                                &
        anlsp, xnlkbsp, nanlsp, lnewanl,                                 &
        mxddim, mxdbnd, mxdasp, mxdgve, mxdscr, mxdnsp)

!   orthogonalization.  Could use SVD...

    allocate(irow(2*neig))

    call grsch_loop_c16(.TRUE., bas, hbas, 2*mtxd, neig, 2*neig, irow,   &
         2*mxddim)

    call grsch_shift_c16(bas, 2*mtxd, neig, 2*neig, irow, ired,          &
         2*mxddim)

    call grsch_shift_c16(hbas, 2*mtxd, neig, 2*neig, irow, ired,         &
         2*mxddim)

    deallocate(irow)

    nbas = 2*neig - ired
    mxdorb = 2*mxdbnd

!   diagonalization

    call diag_psi_h_psi_c16(2*mtxd, neig, eg,                            &
          psi_sp, hpsi_sp, nout, .TRUE.,                                 &
          nbas, bas, hbas,  epsdeg,                                      &
          2*mxddim, mxdbnd, mxdorb)

  enddo

  deallocate(hbas)
  deallocate(bas)

  return

end subroutine diag_jacobi_ritz_spin_c16
